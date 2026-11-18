function [quant, fitData, coordFiles, badCoordTable, cacheEvent, ...
    originalData, preFitBasis] = ...
    ReadPatientCoordCached(coordDir, patientID, cacheOptions, varargin)
%ReadPatientCoordCached Read one patient's LCModel files through one cache.
%
% The cache stores the parsed VDIQuant object, fitData, successfully parsed
% files, and malformed-file diagnostics. With includeRawBasis=true it also
% stores each matching complex .RAW FID and the patient's complex pre-fit
% .basis spectra. Invalidation compares the complete discovered source set
% (ordered paths, names, sizes, and timestamps), the patient/source identity,
% schema version, and parser metadata.
% Malformed sources deliberately remain in that signature so changing or
% replacing one invalidates the patient cache even though it was not parsed.

    p = inputParser;
    p.addParameter('filePrefix', "", @(x) ischar(x) || isstring(x));
    p.addParameter('includeRawBasis', false, @(x) islogical(x) && isscalar(x));
    parse(p, varargin{:});

    filePrefix = string(p.Results.filePrefix);
    includeRawBasis = logical(p.Results.includeRawBasis);
    cacheEvent = MakeEmptyCoordCacheEvent();
    originalData = MakeEmptyOriginalData();
    preFitBasis = MakeEmptyPreFitBasis();

    validationTimer = tic;
    [sourceCoordFiles, coordInfo] = DiscoverCoordFiles(coordDir, filePrefix);

    sourceRawFiles = strings(0, 1);
    rawInfo = struct([]);
    basisFile = "";
    basisInfo = struct([]);
    if includeRawBasis
        [sourceRawFiles, rawInfo, basisFile, basisInfo] = ...
            DiscoverRawBasisFiles(coordDir, patientID, sourceCoordFiles);
    end

    sourceMetadata = table();
    inputSourceMetadata = table();
    parserMetadata = struct();
    rawBasisParserMetadata = struct();
    cacheFile = "";

    if cacheOptions.enabled
        sourceMetadata = BuildCoordSourceMetadata(coordInfo);
        parserMetadata = BuildCoordParserMetadata();
        if includeRawBasis
            inputSourceMetadata = BuildRawBasisSourceMetadata(rawInfo, basisInfo);
            rawBasisParserMetadata = BuildRawBasisParserMetadata();
        end
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
            cacheFile, patientID, coordDir, filePrefix, sourceMetadata, parserMetadata, ...
            includeRawBasis, inputSourceMetadata, rawBasisParserMetadata, ...
            sourceRawFiles, basisFile);
        cacheEvent.validationSeconds = cacheEvent.validationSeconds + ...
            max(toc(validationTimer) - loadSeconds, 0);
        cacheEvent.loadSeconds = loadSeconds;

        if isHit
            quant = cached.quant;
            fitData = cached.fitData;
            coordFiles = string(cached.coordFiles(:));
            badCoordTable = cached.badCoordTable;
            if includeRawBasis
                originalData = cached.originalData;
                preFitBasis = cached.preFitBasis;
            end
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
    if includeRawBasis
        [originalData, preFitBasis] = ReadRawBasisInputs( ...
            patientID, coordFiles, fitData, sourceCoordFiles, sourceRawFiles, basisFile);
    end
    cacheEvent.parseSeconds = toc(parseTimer);
    cacheEvent.parserCallCount = parserCallCount;

    if cacheOptions.enabled
        writeTimer = tic;
        cacheEvent.didAttemptWrite = true;
        coordCache = BuildPatientCoordCacheRecord( ...
            patientID, coordDir, filePrefix, sourceMetadata, parserMetadata, ...
            quant, fitData, coordFiles, badCoordTable, includeRawBasis, ...
            inputSourceMetadata, rawBasisParserMetadata, originalData, preFitBasis);
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


function [originalData, preFitBasis] = ReadRawBasisInputs( ...
    patientID, coordFiles, fitData, sourceCoordFiles, sourceRawFiles, basisFile)

    coordFiles = string(coordFiles(:));
    sourceCoordFiles = string(sourceCoordFiles(:));
    sourceRawFiles = string(sourceRawFiles(:));

    [isMapped, sourceIndex] = ismember( ...
        NormalizeCachePath(coordFiles), NormalizeCachePath(sourceCoordFiles));
    if any(~isMapped)
        error('Could not map every successfully parsed .coord file to its .RAW source.');
    end
    rawFiles = sourceRawFiles(sourceIndex);

    try
        [basisSpectra, basisHeaders, basisPPM, basisGlobalHeader] = ...
            VDILCM.ReadBasis(char(basisFile));
    catch ME
        error('ReadPatientCoordCached:BasisReadFailed', ...
            'Could not read patient basis file %s: %s', basisFile, ME.message);
    end

    nFiles = numel(rawFiles);
    rawFIDCells = cell(nFiles, 1);
    rawHeaderCells = cell(nFiles, 1);
    rawPointCounts = zeros(nFiles, 1);

    for k = 1:nFiles
        try
            [rawFIDCells{k}, rawHeaderCells{k}] = VDILCM.ReadRAW(char(rawFiles(k)));
        catch ME
            error('ReadPatientCoordCached:RawReadFailed', ...
                'Could not read .RAW file corresponding to %s:\n%s\n%s', ...
                coordFiles(k), rawFiles(k), ME.message);
        end
        rawPointCounts(k) = numel(rawFIDCells{k});
    end

    ValidateAcquisitionGrid(coordFiles, rawFiles, rawFIDCells, rawHeaderCells, ...
        fitData, basisFile, basisSpectra, basisPPM, basisGlobalHeader);

    nRawPoints = rawPointCounts(1);
    fids = complex(zeros(nRawPoints, nFiles));
    for k = 1:nFiles
        fids(:, k) = rawFIDCells{k}(:);
    end

    division = zeros(nFiles, 1);
    part = zeros(nFiles, 1);
    for k = 1:nFiles
        [division(k), part(k)] = ParseDivisionPart(coordFiles(k));
    end
    rowIndex = (1:nFiles).';
    fileTable = table(rowIndex, division, part, coordFiles, rawFiles);

    originalData = struct();
    originalData.patientID = string(patientID);
    originalData.fileTable = fileTable;
    originalData.fids = fids;
    originalData.headers = vertcat(rawHeaderCells{:});
    originalData.metadata = struct( ...
        'domain', "time", ...
        'pointDimension', 1, ...
        'fileDimension', 2, ...
        'pointCount', nRawPoints, ...
        'dwellTimeSeconds', double(basisGlobalHeader.badelt), ...
        'spectralWidthHz', 1 / double(basisGlobalHeader.badelt), ...
        'isComplex', ~isreal(fids), ...
        'transformedByLoader', false);

    basisPointCount = numel(basisSpectra{1});
    basisNames = string(arrayfun(@(x) x.metabo, basisHeaders, ...
        'UniformOutput', false));
    preFitBasis = struct();
    preFitBasis.patientID = string(patientID);
    preFitBasis.file = string(basisFile);
    preFitBasis.spectra = basisSpectra;
    preFitBasis.headers = basisHeaders;
    preFitBasis.metaboliteNames = basisNames(:);
    preFitBasis.ppmAxis = basisPPM(:);
    preFitBasis.globalHeader = basisGlobalHeader;
    preFitBasis.metadata = struct( ...
        'domain', "frequency", ...
        'fftShiftedByReadBasis', true, ...
        'pointCount', basisPointCount, ...
        'originalFidPointCount', nRawPoints, ...
        'zeroFillFactor', basisPointCount / nRawPoints, ...
        'dwellTimeSeconds', double(basisGlobalHeader.badelt), ...
        'spectralWidthHz', 1 / double(basisGlobalHeader.badelt), ...
        'ppmAxisSource', "VDILCM.ReadBasis inferred peak reference", ...
        'transformedByLoader', false);
end


function ValidateAcquisitionGrid(coordFiles, rawFiles, rawFIDs, rawHeaders, ...
    fitData, basisFile, basisSpectra, basisPPM, basisGlobalHeader)

    requiredBasisFields = {'badelt', 'hzpppm', 'echot', 'ndatab', 'vectorsize'};
    for k = 1:numel(requiredBasisFields)
        if ~isfield(basisGlobalHeader, requiredBasisFields{k})
            error('Basis file %s lacks required grid field %s.', ...
                basisFile, requiredBasisFields{k});
        end
    end
    if isempty(basisSpectra)
        error('Basis file contains no basis functions: %s', basisFile);
    end

    rawPointCounts = cellfun(@numel, rawFIDs);
    if any(rawPointCounts ~= rawPointCounts(1))
        badIndex = find(rawPointCounts ~= rawPointCounts(1), 1);
        error('ReadPatientCoordCached:InconsistentRawGrid', ...
            ['Inconsistent .RAW point count for %s: expected %d from %s, ', ...
             'found %d.'], rawFiles(badIndex), rawPointCounts(1), ...
            rawFiles(1), rawPointCounts(badIndex));
    end

    nRawPoints = rawPointCounts(1);
    nBasisPoints = cellfun(@numel, basisSpectra);
    if any(nBasisPoints ~= nBasisPoints(1))
        error('Basis functions in %s do not share one point count.', basisFile);
    end
    if numel(basisPPM) ~= nBasisPoints(1)
        error('Basis PPM axis length does not match the basis spectra in %s.', basisFile);
    end
    if double(basisGlobalHeader.vectorsize) ~= nRawPoints
        error('ReadPatientCoordCached:IncompatibleBasisGrid', ...
            ['RAW/basis grid mismatch for patient basis %s: RAW has %d points, ', ...
             'basis vectorsize is %g.'], ...
            basisFile, nRawPoints, double(basisGlobalHeader.vectorsize));
    end
    if double(basisGlobalHeader.ndatab) ~= nBasisPoints(1)
        error('Basis ndatab (%g) does not match parsed basis length (%d) in %s.', ...
            double(basisGlobalHeader.ndatab), nBasisPoints(1), basisFile);
    end

    for k = 1:numel(rawHeaders)
        rawHeader = rawHeaders{k};
        ValidateNearlyEqual(rawHeader.HZPPPM, basisGlobalHeader.hzpppm, ...
            'HZPPPM', rawFiles(k), basisFile);
        ValidateNearlyEqual(rawHeader.ECHOT, basisGlobalHeader.echot, ...
            'echo time', rawFiles(k), basisFile);

        if k <= numel(fitData) && isfield(fitData(k), 'hdr')
            fitHeader = fitData(k).hdr;
            if isfield(fitHeader, 'nunfil') && double(fitHeader.nunfil) ~= nRawPoints
                error('RAW point count does not match coord nunfil for %s.', coordFiles(k));
            end
            if isfield(fitHeader, 'dwellTime')
                ValidateNearlyEqual(fitHeader.dwellTime, basisGlobalHeader.badelt, ...
                    'dwell time', coordFiles(k), basisFile);
            end
            if isfield(fitHeader, 'PPM')
                ValidateNearlyEqual(fitHeader.PPM, basisGlobalHeader.hzpppm, ...
                    'HZPPPM', coordFiles(k), basisFile);
            end
            if isfield(fitHeader, 'TE')
                ValidateNearlyEqual(fitHeader.TE, basisGlobalHeader.echot, ...
                    'echo time', coordFiles(k), basisFile);
            end
        end
    end
end


function ValidateNearlyEqual(valueA, valueB, quantity, fileA, fileB)
    valueA = double(valueA);
    valueB = double(valueB);
    if ~isscalar(valueA) || ~isscalar(valueB) || ...
            ~isfinite(valueA) || ~isfinite(valueB)
        error('ReadPatientCoordCached:IncompatibleAcquisitionGrid', ...
            '%s mismatch between %s (%g) and %s (%g).', ...
            quantity, fileA, valueA, fileB, valueB);
    end
    tolerance = max(1e-10, 1e-6 * max(abs([valueA, valueB])));
    if abs(valueA - valueB) > tolerance
        error('ReadPatientCoordCached:IncompatibleAcquisitionGrid', ...
            '%s mismatch between %s (%g) and %s (%g).', ...
            quantity, fileA, valueA, fileB, valueB);
    end
end


function [division, part] = ParseDivisionPart(filename)
    token = regexp(char(filename), ...
        'Division_(\d+)_(?:part_)?(\d+)\.basis\.coord$', 'tokens', 'once');
    if isempty(token)
        error('Could not parse division/part identity from: %s', filename);
    end
    division = str2double(token{1});
    part = str2double(token{2});
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


function [rawFiles, rawInfo, basisFile, basisInfo] = ...
    DiscoverRawBasisFiles(coordDir, patientID, coordFiles)

    coordDir = string(coordDir);
    coordFiles = string(coordFiles(:));
    rawFiles = regexprep(coordFiles, '\.coord$', '.RAW', 'ignorecase');

    missingRaw = rawFiles(~isfile(rawFiles));
    if ~isempty(missingRaw)
        error('ReadPatientCoordCached:MissingRawFiles', ...
            'Missing %d matching .RAW file(s) for patient %s:\n%s', ...
            numel(missingRaw), string(patientID), strjoin(missingRaw, newline));
    end

    rawInfo = repmat(dir(char(rawFiles(1))), numel(rawFiles), 1);
    for k = 1:numel(rawFiles)
        curInfo = dir(char(rawFiles(k)));
        if numel(curInfo) ~= 1 || curInfo.isdir
            error('Could not resolve exactly one .RAW source: %s', rawFiles(k));
        end
        rawInfo(k) = curInfo;
    end

    basisInfo = dir(fullfile(coordDir, '*.basis'));
    basisInfo = basisInfo(~[basisInfo.isdir]);
    if isempty(basisInfo)
        error('ReadPatientCoordCached:MissingBasisFile', ...
            'No patient-level .basis file found for patient %s in: %s', ...
            string(patientID), coordDir);
    end

    if numel(basisInfo) > 1
        safePatientID = string(matlab.lang.makeValidName(char(string(patientID))));
        expectedPrefix = "For_division_" + safePatientID + "_";
        isPatientBasis = startsWith(string({basisInfo.name}), expectedPrefix, ...
            'IgnoreCase', true);
        if nnz(isPatientBasis) == 1
            basisInfo = basisInfo(isPatientBasis);
        else
            error('ReadPatientCoordCached:AmbiguousBasisFile', ...
                'Expected one patient-level .basis file for patient %s in %s; found:\n%s', ...
                string(patientID), coordDir, ...
                strjoin(string({basisInfo.name}), newline));
        end
    end

    basisFile = string(fullfile(basisInfo.folder, basisInfo.name));
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


function metadata = BuildRawBasisSourceMetadata(rawInfo, basisInfo)
    nRaw = numel(rawInfo);
    fileType = [repmat("raw", nRaw, 1); "basis"];
    orderIndex = (1:(nRaw + 1)).';
    allInfo = [rawInfo(:); basisInfo(:)];
    fullPath = strings(nRaw + 1, 1);
    filename = strings(nRaw + 1, 1);
    bytes = zeros(nRaw + 1, 1);
    modifiedDatenum = zeros(nRaw + 1, 1);

    for k = 1:numel(allInfo)
        fullPath(k) = NormalizeCachePath(fullfile(allInfo(k).folder, allInfo(k).name));
        filename(k) = string(allInfo(k).name);
        bytes(k) = double(allInfo(k).bytes);
        modifiedDatenum(k) = double(allInfo(k).datenum);
    end

    metadata = table(fileType, orderIndex, fullPath, filename, bytes, modifiedDatenum);
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


function metadata = BuildRawBasisParserMetadata()
    parserFile = string(which('VDILCM'));
    if strlength(parserFile) == 0 || ~isfile(parserFile)
        error('Could not locate VDILCM.m for RAW/basis-cache validation.');
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
    sourceMetadata, parserMetadata, includeRawBasis, inputSourceMetadata, ...
    rawBasisParserMetadata, sourceRawFiles, basisFile)

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
        cached, patientID, coordDir, filePrefix, sourceMetadata, parserMetadata, ...
        includeRawBasis, inputSourceMetadata, rawBasisParserMetadata, ...
        sourceRawFiles, basisFile);

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
    cached, patientID, coordDir, filePrefix, sourceMetadata, parserMetadata, ...
    includeRawBasis, inputSourceMetadata, rawBasisParserMetadata, ...
    sourceRawFiles, basisFile)

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

    if includeRawBasis
        extensionFields = ["inputSourceMetadata", "rawBasisParserMetadata", ...
            "originalData", "preFitBasis"];
        for fieldName = extensionFields
            if ~isfield(cached, fieldName)
                reason = "cache does not yet contain requested RAW/basis data";
                return;
            end
        end

        if ~istable(cached.inputSourceMetadata) || ...
                ~isequaln(cached.inputSourceMetadata, inputSourceMetadata)
            reason = "RAW/basis source metadata changed";
            return;
        end
        if ~isstruct(cached.rawBasisParserMetadata) || ...
                ~isequaln(cached.rawBasisParserMetadata, rawBasisParserMetadata)
            reason = "VDILCM parser metadata changed";
            return;
        end

        [extensionIsValid, extensionReason] = ValidateCachedRawBasis( ...
            cached.originalData, cached.preFitBasis, goodFiles, ...
            sourceRawFiles, basisFile);
        if ~extensionIsValid
            reason = extensionReason;
            isCorrupt = true;
            return;
        end
    end

    isValid = true;
end


function [isValid, reason] = ValidateCachedRawBasis( ...
    originalData, preFitBasis, goodFiles, sourceRawFiles, basisFile)

    isValid = false;
    reason = "";
    originalFields = ["fileTable", "fids", "headers", "metadata"];
    basisFields = ["file", "spectra", "headers", "metaboliteNames", ...
        "ppmAxis", "globalHeader", "metadata"];

    if ~isstruct(originalData) || ~isscalar(originalData) || ...
            ~all(isfield(originalData, originalFields))
        reason = "cached originalData record is invalid";
        return;
    end
    if ~isstruct(preFitBasis) || ~isscalar(preFitBasis) || ...
            ~all(isfield(preFitBasis, basisFields))
        reason = "cached preFitBasis record is invalid";
        return;
    end
    if ~istable(originalData.fileTable) || ...
            ~all(ismember(["division", "part", "coordFiles", "rawFiles"], ...
            string(originalData.fileTable.Properties.VariableNames)))
        reason = "cached RAW file table schema is invalid";
        return;
    end

    cachedCoordFiles = NormalizeCachePath(string(originalData.fileTable.coordFiles));
    cachedRawFiles = NormalizeCachePath(string(originalData.fileTable.rawFiles));
    sourceRawFiles = string(sourceRawFiles(:));
    sourceCoordFiles = regexprep(sourceRawFiles, '\.RAW$', '.coord', 'ignorecase');
    [isMapped, sourceIndex] = ismember(goodFiles, NormalizeCachePath(sourceCoordFiles));
    if any(~isMapped)
        reason = "cached RAW mapping does not cover all good coord files";
        return;
    end
    expectedRawFiles = NormalizeCachePath(sourceRawFiles(sourceIndex));

    nFiles = numel(goodFiles);
    if ~isequal(cachedCoordFiles, goodFiles) || ...
            ~isequal(cachedRawFiles, expectedRawFiles) || ...
            ~isnumeric(originalData.fids) || size(originalData.fids, 2) ~= nFiles || ...
            ~isstruct(originalData.headers) || numel(originalData.headers) ~= nFiles
        reason = "cached RAW data/file accounting is inconsistent";
        return;
    end
    if ~iscell(preFitBasis.spectra) || isempty(preFitBasis.spectra) || ...
            ~isstruct(preFitBasis.headers) || ...
            numel(preFitBasis.headers) ~= numel(preFitBasis.spectra) || ...
            ~isnumeric(preFitBasis.ppmAxis) || ...
            numel(preFitBasis.ppmAxis) ~= numel(preFitBasis.spectra{1}) || ...
            ~isstruct(preFitBasis.globalHeader) || ...
            ~strcmp(NormalizeCachePath(preFitBasis.file), NormalizeCachePath(basisFile))
        reason = "cached pre-fit basis data is inconsistent";
        return;
    end

    isValid = true;
end


function cacheRecord = BuildPatientCoordCacheRecord( ...
    patientID, coordDir, filePrefix, sourceMetadata, parserMetadata, ...
    quant, fitData, coordFiles, badCoordTable, includeRawBasis, ...
    inputSourceMetadata, rawBasisParserMetadata, originalData, preFitBasis)

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
    if includeRawBasis
        cacheRecord.inputSourceMetadata = inputSourceMetadata;
        cacheRecord.rawBasisParserMetadata = rawBasisParserMetadata;
        cacheRecord.originalData = originalData;
        cacheRecord.preFitBasis = preFitBasis;
    end
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


function data = MakeEmptyOriginalData()
    data = struct();
    data.patientID = "";
    data.fileTable = table();
    data.fids = complex(zeros(0, 0));
    data.headers = struct([]);
    data.metadata = struct();
end


function basis = MakeEmptyPreFitBasis()
    basis = struct();
    basis.patientID = "";
    basis.file = "";
    basis.spectra = cell(0, 1);
    basis.headers = struct([]);
    basis.metaboliteNames = strings(0, 1);
    basis.ppmAxis = zeros(0, 1);
    basis.globalHeader = struct();
    basis.metadata = struct();
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
