function outputs = MetabCovarianceByPatient(cfg, varargin)
% MetabCovarianceByPatient
%
% Loads LCModel .coord files from patient folders and computes:
%
%   1. A metabolite-by-metabolite covariance matrix for each patient
%   2. A metabolite-by-metabolite correlation matrix for each patient
%   3. The average covariance matrix across patients
%   4. The average correlation matrix across patients
%
% Default behavior:
%   - uses Division_1
%   - uses sig values
%   - loads folders according to cfg.load or cfg.covariance
%
% Main output examples:
%
%   outputs.patientResults(1).covTable
%   outputs.patientResults(1).corrTable
%
%   outputs.patientResultsByID.P11.covTable
%   outputs.patientResultsByID.P11.corrTable
%
%   outputs.group.meanCovTable
%   outputs.group.meanCorrTable

    if nargin < 1 || isempty(cfg)
        cfg = ProjectConfig();

    elseif ischar(cfg) || isstring(cfg)
        varargin = [{cfg}, varargin];
        cfg = ProjectConfig();
    end

    opts = ParseCovarianceOptionsFromConfig(cfg, varargin{:});

    patients = BuildPatientFolderTableForCovariance(opts);

    nPatients = height(patients);
    metabList = string(opts.metabList(:));
    nMetabs = numel(metabList);

    patientResults = struct([]);
    patientResultsByID = struct;

    covStack = nan(nMetabs, nMetabs, nPatients);
    corrStack = nan(nMetabs, nMetabs, nPatients);
    absCorrStack = nan(nMetabs, nMetabs, nPatients);
    nPairStack = nan(nMetabs, nMetabs, nPatients);

    fprintf('\nMetabolite covariance/correlation analysis\n');
    fprintf('Load mode: %s\n', opts.loadMode);
    fprintf('Number of patient folders: %d\n', nPatients);
    fprintf('Division used: Division_%d\n', opts.division);
    fprintf('Value column: %s\n', opts.valueColumn);

    for pIdx = 1:nPatients

        patientID = string(patients.patientID(pIdx));
        coordDir = string(patients.coordDir(pIdx));

        fprintf('\nReading patient %s\n', patientID);
        fprintf('Coord dir: %s\n', coordDir);

        if opts.addPatientPrefixToFilenames
            filePrefix = patientID + "_";
        else
            filePrefix = "";
        end

        [quant, fitData, coordFiles, badCoordTable] = SafeReadCoordFiles( ...
            coordDir, ...
            'filePrefix', filePrefix);

        coordTable = quant.metabTable;

        partTable = BuildDivisionMetabPartTable( ...
            coordTable, ...
            metabList, ...
            opts.valueColumn, ...
            opts.division, ...
            opts.ignoreZeros);

        dataMatrix = table2array(partTable(:, 2:end));

        [covMatrix, corrMatrix, nPairMatrix] = ComputePairwiseCovCorr( ...
            dataMatrix, ...
            opts.minValidPairs);

        absCorrMatrix = abs(corrMatrix);

        covTable = MatrixToMetabTable(covMatrix, metabList);
        corrTable = MatrixToMetabTable(corrMatrix, metabList);
        absCorrTable = MatrixToMetabTable(absCorrMatrix, metabList);
        nPairTable = MatrixToMetabTable(nPairMatrix, metabList);

        patientResults(pIdx).patientID = patientID;
        patientResults(pIdx).coordDir = coordDir;
        patientResults(pIdx).coordFiles = coordFiles;
        patientResults(pIdx).badCoordTable = badCoordTable;
        patientResults(pIdx).coordTable = coordTable;
        patientResults(pIdx).fitData = fitData;

        patientResults(pIdx).division = opts.division;
        patientResults(pIdx).valueColumn = opts.valueColumn;
        patientResults(pIdx).metabList = metabList;

        patientResults(pIdx).partTable = partTable;
        patientResults(pIdx).dataMatrix = dataMatrix;

        patientResults(pIdx).covMatrix = covMatrix;
        patientResults(pIdx).corrMatrix = corrMatrix;
        patientResults(pIdx).absCorrMatrix = absCorrMatrix;
        patientResults(pIdx).nPairMatrix = nPairMatrix;

        patientResults(pIdx).covTable = covTable;
        patientResults(pIdx).corrTable = corrTable;
        patientResults(pIdx).absCorrTable = absCorrTable;
        patientResults(pIdx).nPairTable = nPairTable;

        safeField = matlab.lang.makeValidName(char(patientID));
        patientResultsByID.(safeField) = patientResults(pIdx);

        covStack(:, :, pIdx) = covMatrix;
        corrStack(:, :, pIdx) = corrMatrix;
        absCorrStack(:, :, pIdx) = absCorrMatrix;
        nPairStack(:, :, pIdx) = nPairMatrix;
    end

    meanCovMatrix = mean(covStack, 3, 'omitnan');
    meanCorrMatrix = AverageCorrelationMatrices(corrStack);

    % This is mean(abs(patient correlations)), not abs(mean correlation).
    % It captures correlation magnitude even when signs differ between patients.
    meanAbsCorrMatrix = mean(absCorrStack, 3, 'omitnan');

    nPatientsCovMatrix = sum(~isnan(covStack), 3);
    nPatientsCorrMatrix = sum(~isnan(corrStack), 3);
    nPatientsAbsCorrMatrix = sum(~isnan(absCorrStack), 3);

    outputs = struct;

    outputs.settings = opts;
    outputs.patients = patients;
    outputs.metabList = metabList;

    outputs.patientResults = patientResults;
    outputs.patientResultsByID = patientResultsByID;

    outputs.group.covStack = covStack;
    outputs.group.corrStack = corrStack;
    outputs.group.absCorrStack = absCorrStack;
    outputs.group.nPairStack = nPairStack;

    outputs.group.meanCovMatrix = meanCovMatrix;
    outputs.group.meanCorrMatrix = meanCorrMatrix;
    outputs.group.meanAbsCorrMatrix = meanAbsCorrMatrix;

    outputs.group.nPatientsCovMatrix = nPatientsCovMatrix;
    outputs.group.nPatientsCorrMatrix = nPatientsCorrMatrix;
    outputs.group.nPatientsAbsCorrMatrix = nPatientsAbsCorrMatrix;

    outputs.group.meanCovTable = MatrixToMetabTable(meanCovMatrix, metabList);
    outputs.group.meanCorrTable = MatrixToMetabTable(meanCorrMatrix, metabList);
    outputs.group.meanAbsCorrTable = MatrixToMetabTable(meanAbsCorrMatrix, metabList);

    outputs.group.nPatientsCovTable = MatrixToMetabTable(nPatientsCovMatrix, metabList);
    outputs.group.nPatientsCorrTable = MatrixToMetabTable(nPatientsCorrMatrix, metabList);
    outputs.group.nPatientsAbsCorrTable = MatrixToMetabTable(nPatientsAbsCorrMatrix, metabList);
end


function opts = ParseCovarianceOptionsFromConfig(cfg, varargin)

    opts = struct;

    opts.coordRoot = "";
    opts.loadMode = "allSubfolders";
    opts.selectedPatientID = "";
    opts.selectedPatientIDs = strings(0, 1);
    opts.selectedSubfolder = "";
    opts.selectedSubfolders = strings(0, 1);

    opts.division = 1;
    opts.valueColumn = "sig";
    opts.ignoreZeros = true;
    opts.minValidPairs = 3;

    opts.addPatientPrefixToFilenames = true;

    opts.metabList = DefaultAnalysisMetabList();

    if HasFieldOrProperty(cfg, "paths")
        cfgPaths = GetFieldOrProperty(cfg, "paths");

        if HasFieldOrProperty(cfgPaths, "coordDir")
            opts.coordRoot = string(GetFieldOrProperty(cfgPaths, "coordDir"));
        end
    end

    if HasFieldOrProperty(cfg, "output")
        cfgOutput = GetFieldOrProperty(cfg, "output");

        if HasFieldOrProperty(cfgOutput, "addPatientPrefixToFilenames")
            opts.addPatientPrefixToFilenames = logical(GetFieldOrProperty(cfgOutput, "addPatientPrefixToFilenames"));
        end
    end

    if HasFieldOrProperty(cfg, "load")
        cfgLoad = GetFieldOrProperty(cfg, "load");

        if HasFieldOrProperty(cfgLoad, "mode")
            opts.loadMode = string(GetFieldOrProperty(cfgLoad, "mode"));
        end

        if HasFieldOrProperty(cfgLoad, "coordDir")
            opts.coordRoot = string(GetFieldOrProperty(cfgLoad, "coordDir"));
        end

        if HasFieldOrProperty(cfgLoad, "selectedPatientID")
            opts.selectedPatientID = string(GetFieldOrProperty(cfgLoad, "selectedPatientID"));
        end

        if HasFieldOrProperty(cfgLoad, "selectedPatientIDs")
            opts.selectedPatientIDs = string(GetFieldOrProperty(cfgLoad, "selectedPatientIDs"));
        end

        if HasFieldOrProperty(cfgLoad, "selectedSubfolder")
            opts.selectedSubfolder = string(GetFieldOrProperty(cfgLoad, "selectedSubfolder"));
        end

        if HasFieldOrProperty(cfgLoad, "selectedSubfolders")
            opts.selectedSubfolders = string(GetFieldOrProperty(cfgLoad, "selectedSubfolders"));
        end
    end

    if HasFieldOrProperty(cfg, "covariance")
        cfgCov = GetFieldOrProperty(cfg, "covariance");

        if HasFieldOrProperty(cfgCov, "coordRoot")
            opts.coordRoot = string(GetFieldOrProperty(cfgCov, "coordRoot"));
        end

        if HasFieldOrProperty(cfgCov, "coordDir")
            opts.coordRoot = string(GetFieldOrProperty(cfgCov, "coordDir"));
        end

        if HasFieldOrProperty(cfgCov, "loadMode")
            opts.loadMode = string(GetFieldOrProperty(cfgCov, "loadMode"));
        end

        if HasFieldOrProperty(cfgCov, "selectedPatientID")
            opts.selectedPatientID = string(GetFieldOrProperty(cfgCov, "selectedPatientID"));
        end

        if HasFieldOrProperty(cfgCov, "selectedPatientIDs")
            opts.selectedPatientIDs = string(GetFieldOrProperty(cfgCov, "selectedPatientIDs"));
        end

        if HasFieldOrProperty(cfgCov, "selectedSubfolder")
            opts.selectedSubfolder = string(GetFieldOrProperty(cfgCov, "selectedSubfolder"));
        end

        if HasFieldOrProperty(cfgCov, "selectedSubfolders")
            opts.selectedSubfolders = string(GetFieldOrProperty(cfgCov, "selectedSubfolders"));
        end

        if HasFieldOrProperty(cfgCov, "division")
            opts.division = double(GetFieldOrProperty(cfgCov, "division"));
        end

        if HasFieldOrProperty(cfgCov, "valueColumn")
            opts.valueColumn = string(GetFieldOrProperty(cfgCov, "valueColumn"));
        end

        if HasFieldOrProperty(cfgCov, "ignoreZeros")
            opts.ignoreZeros = logical(GetFieldOrProperty(cfgCov, "ignoreZeros"));
        end

        if HasFieldOrProperty(cfgCov, "minValidPairs")
            opts.minValidPairs = double(GetFieldOrProperty(cfgCov, "minValidPairs"));
        end

        if HasFieldOrProperty(cfgCov, "metabList")
            opts.metabList = string(GetFieldOrProperty(cfgCov, "metabList"));
        end
    end

    opts.metabList = TryGetMetabListFromConfig(cfg, opts.metabList);

    p = inputParser;

    p.addParameter('coordRoot', opts.coordRoot, @(x) ischar(x) || isstring(x));
    p.addParameter('loadMode', opts.loadMode, @(x) ischar(x) || isstring(x));

    p.addParameter('selectedPatientID', opts.selectedPatientID, @(x) ischar(x) || isstring(x));
    p.addParameter('selectedPatientIDs', opts.selectedPatientIDs, @(x) ischar(x) || isstring(x) || iscellstr(x));

    p.addParameter('selectedSubfolder', opts.selectedSubfolder, @(x) ischar(x) || isstring(x));
    p.addParameter('selectedSubfolders', opts.selectedSubfolders, @(x) ischar(x) || isstring(x) || iscellstr(x));

    p.addParameter('division', opts.division, @(x) isnumeric(x) && isscalar(x));
    p.addParameter('valueColumn', opts.valueColumn, @(x) ischar(x) || isstring(x));
    p.addParameter('ignoreZeros', opts.ignoreZeros, @(x) islogical(x) && isscalar(x));
    p.addParameter('minValidPairs', opts.minValidPairs, @(x) isnumeric(x) && isscalar(x));

    p.addParameter('addPatientPrefixToFilenames', opts.addPatientPrefixToFilenames, ...
        @(x) islogical(x) && isscalar(x));

    p.addParameter('metabList', opts.metabList, @(x) ischar(x) || isstring(x) || iscellstr(x));

    parse(p, varargin{:});

    opts.coordRoot = string(p.Results.coordRoot);
    opts.loadMode = string(p.Results.loadMode);

    opts.selectedPatientID = string(p.Results.selectedPatientID);
    opts.selectedPatientIDs = string(p.Results.selectedPatientIDs);

    opts.selectedSubfolder = string(p.Results.selectedSubfolder);
    opts.selectedSubfolders = string(p.Results.selectedSubfolders);

    opts.division = double(p.Results.division);
    opts.valueColumn = string(p.Results.valueColumn);
    opts.ignoreZeros = logical(p.Results.ignoreZeros);
    opts.minValidPairs = double(p.Results.minValidPairs);

    opts.addPatientPrefixToFilenames = logical(p.Results.addPatientPrefixToFilenames);

    opts.metabList = string(p.Results.metabList);
    opts.metabList = opts.metabList(:);

    if strlength(opts.coordRoot) == 0
        error('coordRoot is empty. Set cfg.paths.coordDir or cfg.covariance.coordRoot.');
    end
end


function patients = BuildPatientFolderTableForCovariance(opts)

    coordRoot = string(opts.coordRoot);
    loadMode = lower(string(opts.loadMode));

    if ~isfolder(coordRoot)
        error('Coordinate root folder does not exist: %s', coordRoot);
    end

    patientID = strings(0, 1);
    coordDir = strings(0, 1);

    switch loadMode

        case "allsubfolders"

            info = dir(coordRoot);
            info = info([info.isdir]);

            folderNames = string({info.name});
            keep = folderNames ~= "." & folderNames ~= "..";
            info = info(keep);

            for k = 1:numel(info)

                curFolder = string(fullfile(info(k).folder, info(k).name));
                coordInfo = dir(fullfile(curFolder, "*.coord"));

                if isempty(coordInfo)
                    continue;
                end

                curPatientID = string(matlab.lang.makeValidName(info(k).name));

                patientID(end+1, 1) = curPatientID; %#ok<AGROW>
                coordDir(end+1, 1) = curFolder; %#ok<AGROW>
            end

        case "singlesubfolder"

            selectedSubfolder = string(opts.selectedSubfolder);
            selectedPatientID = string(opts.selectedPatientID);

            if strlength(selectedSubfolder) > 0

                curFolder = selectedSubfolder;

                if ~isfolder(curFolder)
                    error('selectedSubfolder does not exist: %s', curFolder);
                end

                if strlength(selectedPatientID) == 0
                    [~, folderName] = fileparts(char(curFolder));
                    selectedPatientID = string(folderName);
                end

            else

                if strlength(selectedPatientID) == 0
                    error('singleSubfolder mode requires selectedPatientID or selectedSubfolder.');
                end

                selectedPatientID = string(matlab.lang.makeValidName(char(selectedPatientID)));
                curFolder = string(fullfile(coordRoot, selectedPatientID));

                if ~isfolder(curFolder)
                    error('Patient folder does not exist: %s', curFolder);
                end
            end

            coordInfo = dir(fullfile(curFolder, "*.coord"));

            if isempty(coordInfo)
                error('No .coord files found in selected folder: %s', curFolder);
            end

            patientID(end+1, 1) = string(matlab.lang.makeValidName(char(selectedPatientID)));
            coordDir(end+1, 1) = curFolder;

        case "selectedsubfolders"

            selectedPatientIDs = string(opts.selectedPatientIDs);
            selectedSubfolders = string(opts.selectedSubfolders);

            selectedPatientIDs = selectedPatientIDs(:);
            selectedSubfolders = selectedSubfolders(:);

            selectedPatientIDs = selectedPatientIDs(strlength(selectedPatientIDs) > 0);
            selectedSubfolders = selectedSubfolders(strlength(selectedSubfolders) > 0);

            for k = 1:numel(selectedPatientIDs)

                curPatientID = string(matlab.lang.makeValidName(char(selectedPatientIDs(k))));
                curFolder = string(fullfile(coordRoot, curPatientID));

                if ~isfolder(curFolder)
                    warning('Skipping missing patient folder: %s', curFolder);
                    continue;
                end

                coordInfo = dir(fullfile(curFolder, "*.coord"));

                if isempty(coordInfo)
                    warning('Skipping folder with no .coord files: %s', curFolder);
                    continue;
                end

                patientID(end+1, 1) = curPatientID; %#ok<AGROW>
                coordDir(end+1, 1) = curFolder; %#ok<AGROW>
            end

            for k = 1:numel(selectedSubfolders)

                curFolder = string(selectedSubfolders(k));

                if ~isfolder(curFolder)
                    warning('Skipping missing selected folder: %s', curFolder);
                    continue;
                end

                coordInfo = dir(fullfile(curFolder, "*.coord"));

                if isempty(coordInfo)
                    warning('Skipping folder with no .coord files: %s', curFolder);
                    continue;
                end

                [~, folderName] = fileparts(char(curFolder));
                curPatientID = string(matlab.lang.makeValidName(folderName));

                patientID(end+1, 1) = curPatientID; %#ok<AGROW>
                coordDir(end+1, 1) = curFolder; %#ok<AGROW>
            end

        otherwise
            error('Unknown loadMode "%s". Use "allSubfolders", "singleSubfolder", or "selectedSubfolders".', opts.loadMode);
    end

    if isempty(patientID)
        error('No valid patient folders were found.');
    end

    patients = table(patientID, coordDir, ...
        'VariableNames', {'patientID', 'coordDir'});
end


function [quant, fitData, coordFiles, badCoordTable] = SafeReadCoordFiles(coordDir, varargin)

    p = inputParser;
    p.addParameter('filePrefix', "", @(x) ischar(x) || isstring(x));
    parse(p, varargin{:});

    filePrefix = string(p.Results.filePrefix);

    if strlength(filePrefix) > 0
        searchPattern = filePrefix + "*.coord";
    else
        searchPattern = "*.coord";
    end

    coordInfo = dir(fullfile(coordDir, searchPattern));

    if isempty(coordInfo)
        error('No .coord files found in %s with pattern %s.', coordDir, searchPattern);
    end

    fileNames = string({coordInfo.name});
    [~, order] = sort(fileNames);
    coordInfo = coordInfo(order);
    fileNames = fileNames(order);

    expectedFile = false(numel(fileNames), 1);

    for k = 1:numel(fileNames)
        expectedFile(k) = ~isempty(regexp(fileNames(k), ...
            '(^|_)Division_\d+_(?:part_)?\d+\.basis\.coord$', ...
            'once'));
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

        try
            VDIIO.ReadLCMCoord(coordFiles(k));

        catch ME
            isGood(k) = false;
            errorMessages(k) = string(ME.message);

            fprintf('\nBad .coord file found:\n%s\n', coordFiles(k));
            fprintf('Error:\n%s\n\n', ME.message);
        end
    end

    badCoordTable = table( ...
        coordFiles(~isGood), ...
        errorMessages(~isGood), ...
        'VariableNames', {'coordFile', 'errorMessage'});

    if all(~isGood)
        error('All .coord files failed to read in: %s', coordDir);
    end

    coordFiles = coordFiles(isGood);

    warning('Reading only %d good .coord files. Skipped %d bad file(s).', ...
        numel(coordFiles), height(badCoordTable));

    [quant, fitData] = VDIIO.ReadLCMCoord(coordFiles);
end


function partTable = BuildDivisionMetabPartTable(coordTable, metabList, valueColumn, divisionNumber, ignoreZeros)

    metabList = string(metabList(:));
    valueColumn = string(valueColumn);

    requiredCols = ["filename", "name", valueColumn];
    tableCols = string(coordTable.Properties.VariableNames);

    for c = requiredCols
        if ~ismember(c, tableCols)
            error('coordTable is missing required column "%s".', c);
        end
    end

    T = coordTable;
    T.filename = string(T.filename);
    T.name = string(T.name);

    nRows = height(T);
    division = nan(nRows, 1);
    part = nan(nRows, 1);

    for r = 1:nRows

        [~, baseName, ext] = fileparts(T.filename(r));
        curFile = string(baseName) + string(ext);

        tok = regexp(curFile, ...
            '.*Division_(\d+)_(?:part_)?(\d+)\.basis\.coord$', ...
            'tokens', 'once');

        if isempty(tok)
            continue;
        end

        division(r) = str2double(tok{1});
        part(r) = str2double(tok{2});
    end

    T.division = division;
    T.part = part;

    T = T(T.division == divisionNumber, :);
    T = T(ismember(T.name, metabList), :);

    if isempty(T)
        error('No rows found for Division_%d and requested metabolites.', divisionNumber);
    end

    parts = unique(T.part);
    parts = sort(parts(:));

    partTable = table;
    partTable.part = parts;

    for m = 1:numel(metabList)

        metabName = metabList(m);
        colName = matlab.lang.makeValidName(char(metabName));

        values = nan(numel(parts), 1);

        for partIdx = 1:numel(parts)

            curPart = parts(partIdx);
            idx = T.name == metabName & T.part == curPart;

            if sum(idx) == 1
                tmp = T.(char(valueColumn))(idx);
                values(partIdx) = double(tmp(1));

            elseif sum(idx) > 1
                warning('Multiple rows found for %s part %d. Using first.', metabName, curPart);
                tmp = T.(char(valueColumn))(idx);
                values(partIdx) = double(tmp(1));
            end
        end

        if ignoreZeros
            values(values == 0) = NaN;
        end

        partTable.(colName) = values;
    end
end


function [covMatrix, corrMatrix, nPairMatrix] = ComputePairwiseCovCorr(dataMatrix, minValidPairs)

    nMetabs = size(dataMatrix, 2);

    covMatrix = nan(nMetabs, nMetabs);
    corrMatrix = nan(nMetabs, nMetabs);
    nPairMatrix = zeros(nMetabs, nMetabs);

    for a = 1:nMetabs

        x = dataMatrix(:, a);

        for b = 1:nMetabs

            y = dataMatrix(:, b);

            valid = ~isnan(x) & ~isnan(y);
            nValid = sum(valid);

            nPairMatrix(a, b) = nValid;

            if nValid < minValidPairs
                continue;
            end

            xv = x(valid);
            yv = y(valid);

            covXY = sum((xv - mean(xv)) .* (yv - mean(yv))) / (nValid - 1);
            covMatrix(a, b) = covXY;

            sx = std(xv, 0);
            sy = std(yv, 0);

            if sx == 0 || sy == 0
                corrMatrix(a, b) = NaN;
            else
                corrMatrix(a, b) = covXY / (sx * sy);
            end
        end
    end
end


function meanCorrMatrix = AverageCorrelationMatrices(corrStack)

    [nMetabs, ~, ~] = size(corrStack);

    meanCorrMatrix = nan(nMetabs, nMetabs);

    for a = 1:nMetabs

        for b = 1:nMetabs

            vals = squeeze(corrStack(a, b, :));
            vals = vals(~isnan(vals));

            if isempty(vals)
                continue;
            end

            if a == b
                meanCorrMatrix(a, b) = 1;
                continue;
            end

            vals(vals >= 1) = 0.999999;
            vals(vals <= -1) = -0.999999;

            zVals = atanh(vals);
            meanCorrMatrix(a, b) = tanh(mean(zVals));
        end
    end
end


function T = MatrixToMetabTable(M, metabList)

    metabList = string(metabList(:));
    varNames = matlab.lang.makeValidName(cellstr(metabList));

    T = array2table(M, ...
        'VariableNames', varNames, ...
        'RowNames', cellstr(metabList));
end


function metabList = TryGetMetabListFromConfig(cfg, fallbackList)

    metabList = string(fallbackList(:));

    possiblePaths = { ...
        ["metabolites", "analysisMetabList"], ...
        ["metabolites", "metabList"], ...
        ["metab", "analysisMetabList"], ...
        ["analysis", "metabList"], ...
        ["analysisMetabList"], ...
        ["metabList"]};

    for k = 1:numel(possiblePaths)

        path = possiblePaths{k};

        try
            cur = cfg;

            for j = 1:numel(path)
                if HasFieldOrProperty(cur, path(j))
                    cur = GetFieldOrProperty(cur, path(j));
                else
                    cur = [];
                    break;
                end
            end

            if ~isempty(cur)
                metabList = string(cur(:));
                return;
            end

        catch
        end
    end
end


function metabList = DefaultAnalysisMetabList()

    basisMetabList = ["NAA", "NAAG", "Cr", "PCr", "GPC", "PCh", "Glu", "Gln", ...
        "GABA", "GSH", "Tau", "Asc", "Glc", "Ace", "mI", "sI", "Asp", "Lac"];

    sumMetabList = ["GPC+PCh", "NAA+NAAG", "Cr+PCr", "Glu+Gln"];

    metabList = [basisMetabList, sumMetabList];
end


function tf = HasFieldOrProperty(S, name)

    name = char(name);

    if isstruct(S)
        tf = isfield(S, name);
    else
        tf = isprop(S, name);
    end
end


function value = GetFieldOrProperty(S, name)

    name = char(name);

    if isstruct(S)
        value = S.(name);
    else
        value = S.(name);
    end
end