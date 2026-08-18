function [quant, fitData, coordFiles, badCoordTable, cacheEvent] = ...
    ReadPatientCoordCached(coordDir, patientID, cacheOptions, varargin)
%ReadPatientCoordCached Read one patient's LCModel files through one cache.
%
% The cache stores the parsed VDIQuant object, fitData, successfully parsed
% files, and malformed-file diagnostics. Invalidation compares the complete
% discovered source set (ordered paths, names, sizes, and timestamps), the
% patient/source identity, schema version, and VDIIO parser metadata.
% Malformed sources deliberately remain in that signature so changing or
% replacing one invalidates the patient cache even though it was not parsed.

    p = inputParser;
    p.addParameter('filePrefix', "", @(x) ischar(x) || isstring(x));
    parse(p, varargin{:});

    filePrefix = string(p.Results.filePrefix);
    cacheEvent = MakeEmptyCoordCacheEvent();

    validationTimer = tic;
    [sourceCoordFiles, coordInfo] = DiscoverCoordFiles(coordDir, filePrefix);

    sourceMetadata = table();
    parserMetadata = struct();
    cacheFile = "";

    if cacheOptions.enabled
        sourceMetadata = BuildCoordSourceMetadata(coordInfo);
        parserMetadata = BuildCoordParserMetadata();
        cacheFile = BuildPatientCoordCacheFile(cacheOptions.directory, patientID);
    end
    cacheEvent.validationSeconds = toc(validationTimer);

    if ~cacheOptions.enabled
        cacheEvent.category = "disabled";
        cacheEvent.status = "DISABLED";
        cacheEvent.reason = "cache disabled";
    elseif cacheOptions.forceRefresh
        cacheEvent.category = "forced";
        cacheEvent.status = "REFRESH";
        cacheEvent.reason = "forced refresh";
    else
        validationTimer = tic;
        [isHit, cached, category, reason, loadSeconds] = TryLoadPatientCoordCache( ...
            cacheFile, patientID, coordDir, filePrefix, sourceMetadata, parserMetadata);
        cacheEvent.validationSeconds = cacheEvent.validationSeconds + ...
            max(toc(validationTimer) - loadSeconds, 0);
        cacheEvent.loadSeconds = loadSeconds;

        if isHit
            quant = cached.quant;
            fitData = cached.fitData;
            coordFiles = string(cached.coordFiles(:));
            badCoordTable = cached.badCoordTable;
            cacheEvent.category = "hit";
            cacheEvent.status = "HIT";
            cacheEvent.reason = "";
            return;
        end

        cacheEvent.category = category;
        cacheEvent.reason = reason;
        switch category
            case "missing"
                cacheEvent.status = "MISS";
            case "invalidated"
                cacheEvent.status = "INVALIDATED";
            otherwise
                cacheEvent.status = "CORRUPT";
        end
    end

    parseTimer = tic;
    [quant, fitData, coordFiles, badCoordTable, parserCallCount] = ReadCoordSafely( ...
        coordDir, sourceCoordFiles);
    cacheEvent.parseSeconds = toc(parseTimer);
    cacheEvent.parserCallCount = parserCallCount;

    if cacheOptions.enabled
        writeTimer = tic;
        cacheEvent.didAttemptWrite = true;
        coordCache = BuildPatientCoordCacheRecord( ...
            patientID, coordDir, filePrefix, sourceMetadata, parserMetadata, ...
            quant, fitData, coordFiles, badCoordTable);
        cacheEvent.writeSucceeded = WritePatientCoordCacheAtomic(cacheFile, coordCache);
        cacheEvent.writeSeconds = toc(writeTimer);
    end
end


function [quant, fitData, coordFiles, badCoordTable, parserCallCount] = ...
    ReadCoordSafely(coordDir, coordFiles)

    parserCallCount = 1;
    badCoordTable = table(strings(0, 1), strings(0, 1), ...
        'VariableNames', {'coordFile', 'errorMessage'});

    try
        [quant, fitData] = VDIIO.ReadLCMCoord(coordFiles);
        return;
    catch ME
        warning('Batch read failed. Checking .coord files one by one...');
        warning('%s', ME.message);
    end

    isGood = true(numel(coordFiles), 1);
    errorMessages = strings(numel(coordFiles), 1);

    for k = 1:numel(coordFiles)
        parserCallCount = parserCallCount + 1;
        try
            VDIIO.ReadLCMCoord(coordFiles(k));
        catch ME
            isGood(k) = false;
            errorMessages(k) = string(ME.message);
            fprintf('\nBad .coord file found:\n%s\n', coordFiles(k));
            fprintf('Error:\n%s\n\n', ME.message);
        end
    end

    badCoordTable = table(coordFiles(~isGood), errorMessages(~isGood), ...
        'VariableNames', {'coordFile', 'errorMessage'});

    if all(~isGood)
        error('All .coord files failed to read in: %s', coordDir);
    end

    coordFiles = coordFiles(isGood);
    warning('Reading only %d good .coord files. Skipped %d bad file(s).', ...
        numel(coordFiles), height(badCoordTable));

    parserCallCount = parserCallCount + 1;
    [quant, fitData] = VDIIO.ReadLCMCoord(coordFiles);
end


function [coordFiles, coordInfo] = DiscoverCoordFiles(coordDir, filePrefix)
    coordDir = string(coordDir);
    filePrefix = string(filePrefix);

    if strlength(filePrefix) > 0
        searchPattern = filePrefix + "*.coord";
    else
        searchPattern = "*.coord";
    end

    coordInfo = dir(fullfile(coordDir, searchPattern));
    if isempty(coordInfo)
        if strlength(filePrefix) > 0
            error('No .coord files found in: %s with prefix: %s', coordDir, filePrefix);
        else
            error('No .coord files found in: %s', coordDir);
        end
    end

    fileNames = string({coordInfo.name});
    [~, order] = sort(fileNames);
    coordInfo = coordInfo(order);
    fileNames = fileNames(order);

    expectedFile = false(numel(fileNames), 1);
    for k = 1:numel(fileNames)
        expectedFile(k) = ~isempty(regexp(fileNames(k), ...
            '(^|_)Division_\d+_(?:part_)?\d+\.basis\.coord$', 'once'));
    end

    if any(~expectedFile)
        warning('Skipping unexpected .coord files in %s:', coordDir);
        disp(fileNames(~expectedFile)')
    end

    coordInfo = coordInfo(expectedFile);
    if isempty(coordInfo)
        error('No expected Division_*.basis.coord files found in: %s', coordDir);
    end

    coordFiles = strings(numel(coordInfo), 1);
    for k = 1:numel(coordInfo)
        coordFiles(k) = fullfile(coordInfo(k).folder, coordInfo(k).name);
    end
end


function metadata = BuildCoordSourceMetadata(coordInfo)
    nFiles = numel(coordInfo);
    orderIndex = (1:nFiles).';
    fullPath = strings(nFiles, 1);
    filename = strings(nFiles, 1);
    bytes = zeros(nFiles, 1);
    modifiedDatenum = zeros(nFiles, 1);

    for k = 1:nFiles
        fullPath(k) = NormalizeCachePath(fullfile(coordInfo(k).folder, coordInfo(k).name));
        filename(k) = string(coordInfo(k).name);
        bytes(k) = double(coordInfo(k).bytes);
        modifiedDatenum(k) = double(coordInfo(k).datenum);
    end

    metadata = table(orderIndex, fullPath, filename, bytes, modifiedDatenum);
end


function metadata = BuildCoordParserMetadata()
    parserFile = string(which('VDIIO'));
    if strlength(parserFile) == 0 || ~isfile(parserFile)
        error('Could not locate VDIIO.m for coord-cache validation.');
    end

    info = dir(parserFile);
    metadata = struct();
    metadata.fullPath = NormalizeCachePath(parserFile);
    metadata.bytes = double(info.bytes);
    metadata.modifiedDatenum = double(info.datenum);
end


function cacheFile = BuildPatientCoordCacheFile(cacheDir, patientID)
    safePatientID = string(matlab.lang.makeValidName( ...
        char(string(patientID)), 'ReplacementStyle', 'hex'));
    cacheFile = string(fullfile(cacheDir, safePatientID + "_coord_cache.mat"));
end


function [isHit, cached, category, reason, loadSeconds] = ...
    TryLoadPatientCoordCache(cacheFile, patientID, coordDir, filePrefix, ...
    sourceMetadata, parserMetadata)

    isHit = false;
    cached = struct();
    category = "missing";
    reason = "no cache file";
    loadSeconds = 0;

    if ~isfile(cacheFile)
        return;
    end

    loadTimer = tic;
    try
        loaded = load(cacheFile, 'coordCache');
    catch ME
        loadSeconds = toc(loadTimer);
        category = "corrupt";
        reason = "cache could not be loaded: " + string(ME.message);
        warning('Ignoring invalid coord cache %s: %s', cacheFile, ME.message);
        return;
    end
    loadSeconds = toc(loadTimer);

    if ~isfield(loaded, 'coordCache') || ~isstruct(loaded.coordCache) || ...
            ~isscalar(loaded.coordCache)
        category = "corrupt";
        reason = "required coordCache record is missing";
        warning('Ignoring invalid coord cache %s: required coordCache record is missing.', cacheFile);
        return;
    end

    cached = loaded.coordCache;
    [isValid, invalidReason, isCorrupt] = ValidatePatientCoordCache( ...
        cached, patientID, coordDir, filePrefix, sourceMetadata, parserMetadata);

    if ~isValid
        if isCorrupt
            category = "corrupt";
            warning('Ignoring invalid coord cache %s: %s', cacheFile, invalidReason);
        else
            category = "invalidated";
        end
        reason = invalidReason;
        cached = struct();
        return;
    end

    isHit = true;
    category = "hit";
    reason = "";
end


function [isValid, reason, isCorrupt] = ValidatePatientCoordCache( ...
    cached, patientID, coordDir, filePrefix, sourceMetadata, parserMetadata)

    isValid = false;
    reason = "";
    isCorrupt = false;
    requiredFields = ["schemaVersion", "patientID", "coordDir", ...
        "filePrefix", "sourceFileCount", "sourceMetadata", "parserMetadata", ...
        "quant", "fitData", "coordFiles", "badCoordTable"];

    for fieldName = requiredFields
        if ~isfield(cached, fieldName)
            reason = "required field is missing: " + fieldName;
            isCorrupt = true;
            return;
        end
    end

    if ~isequal(cached.schemaVersion, CoordCacheSchemaVersion())
        reason = "cache schema version changed";
        return;
    end

    if ~strcmp(string(cached.patientID), string(patientID)) || ...
            ~strcmp(NormalizeCachePath(cached.coordDir), NormalizeCachePath(coordDir)) || ...
            ~strcmp(string(cached.filePrefix), string(filePrefix))
        reason = "patient identity or source directory changed";
        return;
    end

    if ~isscalar(cached.sourceFileCount) || ...
            double(cached.sourceFileCount) ~= height(sourceMetadata)
        reason = "source file count changed";
        return;
    end

    if ~istable(cached.sourceMetadata) || ~isequaln(cached.sourceMetadata, sourceMetadata)
        reason = "source metadata changed";
        return;
    end

    if ~isstruct(cached.parserMetadata) || ~isequaln(cached.parserMetadata, parserMetadata)
        reason = "VDIIO parser metadata changed";
        return;
    end

    if ~isa(cached.quant, 'VDIQuant') || ~isprop(cached.quant, 'metabTable') || ...
            ~istable(cached.quant.metabTable)
        reason = "cached quant object is invalid";
        isCorrupt = true;
        return;
    end

    if ~isstruct(cached.fitData) || ~istable(cached.badCoordTable)
        reason = "cached fitData or badCoordTable has an invalid type";
        isCorrupt = true;
        return;
    end

    badVars = string(cached.badCoordTable.Properties.VariableNames);
    if ~all(ismember(["coordFile", "errorMessage"], badVars))
        reason = "cached badCoordTable schema is invalid";
        isCorrupt = true;
        return;
    end

    goodFiles = NormalizeCachePath(string(cached.coordFiles(:)));
    badFiles = NormalizeCachePath(string(cached.badCoordTable.coordFile(:)));
    expectedGoodFiles = sourceMetadata.fullPath(~ismember(sourceMetadata.fullPath, badFiles));

    if ~isequal(goodFiles, expectedGoodFiles) || ...
            numel(goodFiles) + numel(badFiles) ~= height(sourceMetadata) || ...
            numel(cached.fitData) ~= numel(goodFiles)
        reason = "cached good/bad file accounting is inconsistent";
        isCorrupt = true;
        return;
    end

    isValid = true;
end


function cacheRecord = BuildPatientCoordCacheRecord( ...
    patientID, coordDir, filePrefix, sourceMetadata, parserMetadata, ...
    quant, fitData, coordFiles, badCoordTable)

    cacheRecord = struct();
    cacheRecord.schemaVersion = CoordCacheSchemaVersion();
    cacheRecord.patientID = string(patientID);
    cacheRecord.coordDir = NormalizeCachePath(coordDir);
    cacheRecord.filePrefix = string(filePrefix);
    cacheRecord.sourceFileCount = height(sourceMetadata);
    cacheRecord.sourceMetadata = sourceMetadata;
    cacheRecord.parserMetadata = parserMetadata;
    cacheRecord.coordFiles = string(coordFiles(:));
    cacheRecord.badCoordTable = badCoordTable;
    cacheRecord.quant = quant;
    cacheRecord.fitData = fitData;
    cacheRecord.createdAtUTC = char(datetime('now', 'TimeZone', 'UTC', ...
        'Format', 'yyyy-MM-dd''T''HH:mm:ss.SSSXXX'));
end


function wasWritten = WritePatientCoordCacheAtomic(cacheFile, coordCache)
    wasWritten = false;
    cacheFile = string(cacheFile);
    cacheDir = string(fileparts(cacheFile));

    try
        if ~isfolder(cacheDir)
            mkdir(cacheDir);
        end

        tempCacheFile = string(tempname(cacheDir)) + ".mat";
        cleanupObj = onCleanup(@() DeleteCoordCacheTempFile(tempCacheFile));
        save(char(tempCacheFile), 'coordCache', '-v7.3');

        [moveSucceeded, moveMessage] = movefile(tempCacheFile, cacheFile, 'f');
        if ~moveSucceeded
            error('Could not finalize cache file: %s', moveMessage);
        end
        wasWritten = true;
    catch ME
        warning('Could not write coord cache %s. Continuing without cache: %s', ...
            cacheFile, ME.message);
    end
end


function DeleteCoordCacheTempFile(tempCacheFile)
    if strlength(string(tempCacheFile)) > 0 && isfile(tempCacheFile)
        delete(tempCacheFile);
    end
end


function pathOut = NormalizeCachePath(pathIn)
    pathOut = string(pathIn);
    pathOut = pathOut(:);
    if ispc
        pathOut = replace(pathOut, '/', '\');
        pathOut = lower(pathOut);
    end
end


function version = CoordCacheSchemaVersion()
    version = 1;
end


function event = MakeEmptyCoordCacheEvent()
    event = struct();
    event.category = "";
    event.status = "";
    event.reason = "";
    event.validationSeconds = 0;
    event.loadSeconds = 0;
    event.parseSeconds = 0;
    event.writeSeconds = 0;
    event.coordStageSeconds = 0;
    event.parserCallCount = 0;
    event.didAttemptWrite = false;
    event.writeSucceeded = false;
end
