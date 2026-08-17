function processedPatients = Splice_data_multi_patient_input(cfg, varargin)
% Splice_data_multi_patient_input
%
% Multi-patient LCModel fitting pipeline.
%
% Input selection can come from ProjectConfig or from name-value options.
%
% Directory mode:
%   Splice_data_multi_patient_input('inputMode', "directory")
%
% Single-file mode:
%   Splice_data_multi_patient_input('inputMode', "singleFile")
%
% Useful overrides:
%   'dataDir'        directory to search in directory mode
%   'twixFile'       file to use in singleFile mode
%   'filePattern'    usually "*.dat"
%   'recursive'      true/false for directory search
%   'mustContain'    string(s) that filenames must contain
%   'mustNotContain' string(s) that filenames must not contain
    
if nargin < 1 || isempty(cfg)
    cfg = ProjectConfig();

elseif ischar(cfg) || isstring(cfg)
    varargin = [{cfg}, varargin];
    cfg = ProjectConfig();
end

runOptions = ParseRunOptionsFromConfig(cfg, varargin{:});
patients = BuildPatientTableFromConfig(cfg, runOptions);
    
usePatientSubfolders = runOptions.usePatientSubfolders;
addPatientPrefixToFilenames = runOptions.addPatientPrefixToFilenames;
deleteOldLCModelFiles = runOptions.deleteOldLCModelFiles;

patients.resolvedCoordDir = strings(height(patients), 1);
patients.filePrefix = strings(height(patients), 1);

fprintf('\nInput mode: %s\n', runOptions.inputMode);
fprintf('Number of patients/files to process: %d\n', height(patients));

%% Run all patients

for pIdx = 1:height(patients)

    patientID = patients.patientID(pIdx);
    twixFile = patients.twixFile(pIdx);
    baseCoordDir = patients.coordDir(pIdx);

    safePatientID = string(matlab.lang.makeValidName(char(patientID)));
    coordDir = ResolvePatientCoordDir(baseCoordDir, patientID, usePatientSubfolders);

    filePrefix = "";
    if addPatientPrefixToFilenames
        filePrefix = safePatientID + "_";
    end

    patients.resolvedCoordDir(pIdx) = coordDir;
    patients.filePrefix(pIdx) = filePrefix;

    fprintf('\nProcessing patient %s\n', patientID);
    fprintf('Twix file: %s\n', twixFile);
    fprintf('Base output dir: %s\n', baseCoordDir);
    fprintf('Resolved patient output dir: %s\n', coordDir);
    fprintf('Filename prefix: %s\n', filePrefix);

    if strlength(twixFile) == 0
        error('Patient %s is missing a twixFile path.', patientID);
    end

    if strlength(baseCoordDir) == 0
        error('Patient %s is missing a coordDir path.', patientID);
    end

    EnsureFolderExists(coordDir);

    if deleteOldLCModelFiles
        CleanLCModelOutputDir(coordDir);
    end

    %% Load and preprocess

    img = VDIIO.LoadTwix(twixFile, 'isICEChop', cfg.preprocessing.isICEChop);

    img.AddCoils;
    img.ChopPts("numPts", cfg.preprocessing.numPts);
    img.FT;
    img.Phase; % This will automatically phase the data

    %% Create and export basis function for this patient

    [basisFunc, TE] = createBasis(img, cfg.metabolites.basis);

    safeImageName = string(matlab.lang.makeValidName(char(erase(string(img.name), " "))));

    basisName = fullfile(coordDir, ...
        "For_division_" + safePatientID + "_" + safeImageName + ".basis");

    basisFunc.ExportBasisToLCModel(basisName, 'TE', TE);

    %% Splice data to subsets and fit LCModel

    for nAvg = cfg.subsets.setSizes
        FitSubset(img, nAvg, basisName, coordDir, filePrefix);
    end
end

processedPatients = patients;
end

%% Functions

function [myBasis, TE] = createBasis(inputImg, metList)
% createBasis
%
% Creates a simple basis set using the metabolite list in ProjectConfig.
% Note: this currently uses an ideal PRESS sequence, as in your original
% script. If your acquisition is sLASER, discuss with your mentor whether
% the basis simulation should be changed.

    spins = SpinsJ('B0', inputImg.B0);
    spins.AddMetab(metList);

    TE = VDIIO.GetTwixHeaderReport(inputImg.metadata.hdr).TE;

    seq = Sequence.GetIdealSequence("PRESS", "TE", TE);
    [~, TT] = seq.Apply(spins, "isVerbose", true);

    myBasis = VDIBasis(TT, ...
        'B0', inputImg.B0, ...
        'numAcqPts', inputImg.numSpecPts, ...
        'dwellTime', inputImg.dwellTime);
end

function imgBlocks = FitSubset(inputImg, nAvg, basisName, outputDir, filePrefix)
% FitSubset
%
% For a requested division size nAvg:
%   1. Slice consecutive blocks of nAvg sets.
%   2. Average each block over the set dimension.
%   3. Fit each averaged block with LCModel.
%
% Output files are named:
%   <filePrefix>Division_<nAvg>_part_<k>.basis.coord
%
% With the default multi-patient settings, this becomes:
%   P01_Division_6_part_1.basis.coord
% inside:
%   <base LCMFit folder>/P01

    if nargin < 5
        filePrefix = "";
    end

    filePrefix = string(filePrefix);

    setString = "set";
    dimTypes = string(inputImg.GetDimType);
    idxSetDim = find(strcmpi(dimTypes, setString), 1);

    if isempty(idxSetDim)
        error('Could not find dimension "%s". Available dimensions are: %s', ...
            setString, strjoin(dimTypes, ", "));
    end

    totalSets = size(inputImg.data, idxSetDim);

    if mod(totalSets, nAvg) ~= 0
        error(['nAvg must evenly divide the total number of sets. ' ...
            'totalSets = %d, nAvg = %d.'], totalSets, nAvg);
    end

    nBlocks = totalSets / nAvg;
    imgBlocks = cell(nBlocks, 1);

    for k = 1:nBlocks

        firstSet = (k - 1) * nAvg + 1;
        lastSet  = k * nAvg;

        sliceArgs = repmat({':'}, 1, numel(dimTypes));
        sliceArgs{idxSetDim} = sprintf('%d:%d', firstSet, lastSet);

        imgBlock = inputImg.Slice(sliceArgs{:});
        imgBlock = imgBlock.Average("dim", setString);

        fitFileName = filePrefix + "Division_" + string(nAvg) + ...
            "_part_" + string(k) + ".basis";

        imgBlock.FitLCModel(basisName, ...
            'isVerbose', true, ...
            'outputFilename', fitFileName, ...
            'outputDir', outputDir);

        imgBlocks{k} = imgBlock;
    end
end

function CleanLCModelOutputDir(coordDir)

    patterns = ["*.coord", "*.control", "*.print", "*.RAW", "*.ps"];

    for i = 1:numel(patterns)
        delete(fullfile(coordDir, patterns(i)));
    end
end

function EnsureFolderExists(folderPath)

    if ~isfolder(folderPath)
        mkdir(folderPath);
    end
end

function coordDir = ResolvePatientCoordDir(baseCoordDir, patientID, usePatientSubfolders)
% ResolvePatientCoordDir
%
% If usePatientSubfolders is true, returns:
%   fullfile(baseCoordDir, safePatientID)
%
% But if baseCoordDir already ends in the patient ID, it is left unchanged.
% This prevents accidentally creating LCMFit/P01/P01.

    baseCoordDir = string(baseCoordDir);
    patientID = string(patientID);
    safePatientID = string(matlab.lang.makeValidName(char(patientID)));

    if ~usePatientSubfolders
        coordDir = baseCoordDir;
        return;
    end

    [~, lastFolderName] = fileparts(char(baseCoordDir));
    lastFolderName = string(lastFolderName);

    if strcmpi(lastFolderName, patientID) || strcmpi(lastFolderName, safePatientID)
        coordDir = baseCoordDir;
    else
        coordDir = string(fullfile(baseCoordDir, safePatientID));
    end
end

function runOptions = ParseRunOptionsFromConfig(cfg, varargin)
% ParseRunOptionsFromConfig
%
% Reads input/output options from ProjectConfig and lets name-value inputs
% override them.

    runOptions = struct;

    runOptions.inputMode = "directory";
    runOptions.dataDir = "";
    runOptions.twixFile = "";
    runOptions.singleFileName = "";
    runOptions.singlePatientID = "";
    
    runOptions.filePattern = "*.dat";
    runOptions.recursive = false;
    runOptions.mustContain = strings(0, 1);
    runOptions.mustNotContain = strings(0, 1);
    runOptions.patientIDMode = "sequential";
    
    runOptions.coordDir = "";
    
    runOptions.usePatientSubfolders = true;
    runOptions.addPatientPrefixToFilenames = true;
    runOptions.deleteOldLCModelFiles = true;

    if HasFieldOrProperty(cfg, "paths")
        cfgPaths = GetFieldOrProperty(cfg, "paths");
        if HasFieldOrProperty(cfgPaths, "dataDir")
            runOptions.dataDir = string(GetFieldOrProperty(cfgPaths, "dataDir"));
        end
        if HasFieldOrProperty(cfgPaths, "coordDir")
            runOptions.coordDir = string(GetFieldOrProperty(cfgPaths, "coordDir"));
        end
    end

    if HasFieldOrProperty(cfg, "files")
        cfgFiles = GetFieldOrProperty(cfg, "files");
        if HasFieldOrProperty(cfgFiles, "twixFile")
            runOptions.twixFile = string(GetFieldOrProperty(cfgFiles, "twixFile"));
        end
    end

    if HasFieldOrProperty(cfg, "input")
        cfgInput = GetFieldOrProperty(cfg, "input");

        if HasFieldOrProperty(cfgInput, "singleFileName")
            runOptions.singleFileName = string(GetFieldOrProperty(cfgInput, "singleFileName"));
        end

        if HasFieldOrProperty(cfgInput, "mode")
            runOptions.inputMode = string(GetFieldOrProperty(cfgInput, "mode"));
        end
        if HasFieldOrProperty(cfgInput, "directory")
            runOptions.dataDir = string(GetFieldOrProperty(cfgInput, "directory"));
        end
        if HasFieldOrProperty(cfgInput, "singleFile")
            runOptions.twixFile = string(GetFieldOrProperty(cfgInput, "singleFile"));
        end
        if HasFieldOrProperty(cfgInput, "filePattern")
            runOptions.filePattern = string(GetFieldOrProperty(cfgInput, "filePattern"));
        end
        if HasFieldOrProperty(cfgInput, "recursive")
            runOptions.recursive = logical(GetFieldOrProperty(cfgInput, "recursive"));
        end
        if HasFieldOrProperty(cfgInput, "filenameMustContain")
            runOptions.mustContain = string(GetFieldOrProperty(cfgInput, "filenameMustContain"));
        end
        if HasFieldOrProperty(cfgInput, "filenameMustNotContain")
            runOptions.mustNotContain = string(GetFieldOrProperty(cfgInput, "filenameMustNotContain"));
        end
        if HasFieldOrProperty(cfgInput, "patientIDMode")
            runOptions.patientIDMode = string(GetFieldOrProperty(cfgInput, "patientIDMode"));
        end
        if HasFieldOrProperty(cfgInput, "singlePatientID")
            runOptions.singlePatientID = string(GetFieldOrProperty(cfgInput, "singlePatientID"));
        end
    end

    if HasFieldOrProperty(cfg, "output")
        cfgOutput = GetFieldOrProperty(cfg, "output");
        if HasFieldOrProperty(cfgOutput, "usePatientSubfolders")
            runOptions.usePatientSubfolders = logical(GetFieldOrProperty(cfgOutput, "usePatientSubfolders"));
        end
        if HasFieldOrProperty(cfgOutput, "addPatientPrefixToFilenames")
            runOptions.addPatientPrefixToFilenames = logical(GetFieldOrProperty(cfgOutput, "addPatientPrefixToFilenames"));
        end
        if HasFieldOrProperty(cfgOutput, "deleteOldLCModelFiles")
            runOptions.deleteOldLCModelFiles = logical(GetFieldOrProperty(cfgOutput, "deleteOldLCModelFiles"));
        end
    end

    p = inputParser;
    p.KeepUnmatched = true;
    p.addParameter('inputMode', runOptions.inputMode, @(x) ischar(x) || isstring(x));
    p.addParameter('dataDir', runOptions.dataDir, @(x) ischar(x) || isstring(x));
    p.addParameter('twixFile', runOptions.twixFile, @(x) ischar(x) || isstring(x));
    p.addParameter('filePattern', runOptions.filePattern, @(x) ischar(x) || isstring(x));
    p.addParameter('recursive', runOptions.recursive, @(x) islogical(x) && isscalar(x));
    p.addParameter('mustContain', runOptions.mustContain, @(x) ischar(x) || isstring(x) || iscellstr(x));
    p.addParameter('mustNotContain', runOptions.mustNotContain, @(x) ischar(x) || isstring(x) || iscellstr(x));
    p.addParameter('patientIDMode', runOptions.patientIDMode, @(x) ischar(x) || isstring(x));
    p.addParameter('singlePatientID', runOptions.singlePatientID, @(x) ischar(x) || isstring(x));
    p.addParameter('coordDir', runOptions.coordDir, @(x) ischar(x) || isstring(x));
    p.addParameter('usePatientSubfolders', runOptions.usePatientSubfolders, @(x) islogical(x) && isscalar(x));
    p.addParameter('addPatientPrefixToFilenames', runOptions.addPatientPrefixToFilenames, @(x) islogical(x) && isscalar(x));
    p.addParameter('deleteOldLCModelFiles', runOptions.deleteOldLCModelFiles, @(x) islogical(x) && isscalar(x));
    p.addParameter('singleFileName', runOptions.singleFileName, @(x) ischar(x) || isstring(x));
    parse(p, varargin{:});

    runOptions.inputMode = string(p.Results.inputMode);
    runOptions.dataDir = string(p.Results.dataDir);
    runOptions.twixFile = string(p.Results.twixFile);
    runOptions.filePattern = string(p.Results.filePattern);
    runOptions.recursive = logical(p.Results.recursive);
    runOptions.mustContain = string(p.Results.mustContain);
    runOptions.mustNotContain = string(p.Results.mustNotContain);
    runOptions.patientIDMode = string(p.Results.patientIDMode);
    runOptions.singlePatientID = string(p.Results.singlePatientID);
    runOptions.coordDir = string(p.Results.coordDir);
    runOptions.usePatientSubfolders = logical(p.Results.usePatientSubfolders);
    runOptions.addPatientPrefixToFilenames = logical(p.Results.addPatientPrefixToFilenames);
    runOptions.deleteOldLCModelFiles = logical(p.Results.deleteOldLCModelFiles);
    runOptions.singleFileName = string(p.Results.singleFileName);
    
    if strlength(runOptions.twixFile) == 0 && strlength(runOptions.singleFileName) > 0
        runOptions.twixFile = string(fullfile(runOptions.dataDir, runOptions.singleFileName));
    end

end

function patients = BuildPatientTableFromConfig(cfg, runOptions)
% BuildPatientTableFromConfig
%
% Returns a standardized table with:
%   patientID | twixFile | coordDir
%
% Priority:
%   1. cfg.patients, if provided manually.
%   2. cfg.input.mode / name-value inputMode.
%   3. old fallback cfg.files.twixFile.

    if HasFieldOrProperty(cfg, "patients")
        rawPatients = GetFieldOrProperty(cfg, "patients");
        patients = NormalizePatientTable(rawPatients, runOptions.coordDir);
        return;
    end

    inputMode = lower(string(runOptions.inputMode));

    switch inputMode
        case "directory"
            twixFiles = FindRelevantTwixFiles( ...
                runOptions.dataDir, ...
                runOptions.filePattern, ...
                runOptions.recursive, ...
                runOptions.mustContain, ...
                runOptions.mustNotContain);

            patientID = MakePatientIDsFromFiles(twixFiles, runOptions.patientIDMode);
            coordDir = repmat(runOptions.coordDir, numel(twixFiles), 1);
            patients = table(patientID, twixFiles, coordDir, ...
                'VariableNames', {'patientID', 'twixFile', 'coordDir'});

        case {"singlefile", "single_file", "file"}
            patientID = string(runOptions.singlePatientID);
            twixFile = string(runOptions.twixFile);
            coordDir = string(runOptions.coordDir);

            if strlength(twixFile) == 0
                error('singleFile mode was requested, but no twixFile was provided.');
            end

            patients = table(patientID, twixFile, coordDir, ...
                'VariableNames', {'patientID', 'twixFile', 'coordDir'});

        otherwise
            error('Unknown inputMode "%s". Use "directory" or "singleFile".', runOptions.inputMode);
    end
end

function twixFiles = FindRelevantTwixFiles(dataDir, filePattern, recursive, mustContain, mustNotContain)
% FindRelevantTwixFiles
%
% Searches a directory for Twix files and filters by filename strings.

    dataDir = string(dataDir);
    filePattern = string(filePattern);
    mustContain = CleanStringFilter(mustContain);
    mustNotContain = CleanStringFilter(mustNotContain);

    if strlength(dataDir) == 0 || ~isfolder(dataDir)
        error('Input directory does not exist: %s', dataDir);
    end

    if recursive
        info = dir(fullfile(dataDir, "**", filePattern));
    else
        info = dir(fullfile(dataDir, filePattern));
    end

    info = info(~[info.isdir]);

    if isempty(info)
        error('No files found in %s matching pattern %s.', dataDir, filePattern);
    end

    fullPaths = strings(numel(info), 1);
    names = strings(numel(info), 1);

    for k = 1:numel(info)
        fullPaths(k) = string(fullfile(info(k).folder, info(k).name));
        names(k) = string(info(k).name);
    end

    keep = true(numel(fullPaths), 1);
    lowerNames = lower(names);

    for i = 1:numel(mustContain)
        keep = keep & contains(lowerNames, lower(mustContain(i)));
    end

    for i = 1:numel(mustNotContain)
        keep = keep & ~contains(lowerNames, lower(mustNotContain(i)));
    end

    twixFiles = fullPaths(keep);
    twixFiles = sort(twixFiles(:));

    if isempty(twixFiles)
        error(['Files were found in %s, but none passed the relevance filters. ' ...
            'Check filenameMustContain / filenameMustNotContain.'], dataDir);
    end
end

function filters = CleanStringFilter(filters)

    filters = string(filters(:));
    filters = filters(~ismissing(filters));
    filters = filters(strlength(filters) > 0);
end

function patientIDs = MakePatientIDsFromFiles(twixFiles, patientIDMode)

    twixFiles = string(twixFiles(:));
    patientIDMode = lower(string(patientIDMode));
    nFiles = numel(twixFiles);

    patientIDs = strings(nFiles, 1);

    switch patientIDMode
        case "filename"
            for k = 1:nFiles
                [~, baseName] = fileparts(char(twixFiles(k)));
                patientIDs(k) = string(matlab.lang.makeValidName(baseName));
            end

        case "midfid"
            for k = 1:nFiles
                [~, baseName] = fileparts(char(twixFiles(k)));
                tok = regexp(baseName, '(MID\d+).*?(FID\d+)', 'tokens', 'once');
                if isempty(tok)
                    patientIDs(k) = "P" + compose("%02d", k);
                else
                    patientIDs(k) = string(tok{1}) + "_" + string(tok{2});
                end
            end

        otherwise
            patientIDs = "P" + compose("%02d", (1:nFiles)');
    end
end

function patients = NormalizePatientTable(rawPatients, defaultCoordDir)
% NormalizePatientTable
%
% Converts cfg.patients to a table with standard variable names:
%   patientID | twixFile | coordDir

    if istable(rawPatients)

        patientID = ReadFirstAvailableTableVar(rawPatients, ...
            ["patientID", "PatientID", "id", "ID", "name", "Name"]);

        twixFile = ReadFirstAvailableTableVar(rawPatients, ...
            ["twixFile", "TwixFile", "twix", "dataFile", "filename"]);

        coordDir = ReadFirstAvailableTableVar(rawPatients, ...
            ["coordDir", "CoordDir", "lcmDir", "LCMDir", "outputDir"]);

    elseif isstruct(rawPatients)

        patientID = ReadFirstAvailableStructField(rawPatients, ...
            ["patientID", "PatientID", "id", "ID", "name", "Name"]);

        twixFile = ReadFirstAvailableStructField(rawPatients, ...
            ["twixFile", "TwixFile", "twix", "dataFile", "filename"]);

        coordDir = ReadFirstAvailableStructField(rawPatients, ...
            ["coordDir", "CoordDir", "lcmDir", "LCMDir", "outputDir"]);

    else
        error('cfg.patients must be either a table or a struct array.');
    end

    patientID = string(patientID(:));
    twixFile = string(twixFile(:));
    coordDir = string(coordDir(:));

    if isempty(twixFile)
        error('Could not find twixFile/dataFile/filename in cfg.patients.');
    end

    if isempty(patientID)
        patientID = "P" + compose("%02d", (1:numel(twixFile))');
    end

    if isempty(coordDir)
        coordDir = repmat(string(defaultCoordDir), numel(twixFile), 1);
    end

    if numel(patientID) ~= numel(twixFile)
        error('patientID and twixFile must have the same number of rows.');
    end

    if numel(coordDir) ~= numel(twixFile)
        error('coordDir and twixFile must have the same number of rows.');
    end

    patients = table(patientID, twixFile, coordDir, ...
        'VariableNames', {'patientID', 'twixFile', 'coordDir'});
end

function value = ReadFirstAvailableTableVar(T, candidateNames)

    candidateNames = string(candidateNames);
    tableVars = string(T.Properties.VariableNames);

    value = strings(0, 1);

    for i = 1:numel(candidateNames)
        idx = find(strcmpi(tableVars, candidateNames(i)), 1);
        if ~isempty(idx)
            value = T.(tableVars(idx));
            return;
        end
    end
end

function value = ReadFirstAvailableStructField(S, candidateNames)

    candidateNames = string(candidateNames);
    structFields = string(fieldnames(S));

    value = strings(0, 1);

    for i = 1:numel(candidateNames)
        idx = find(strcmpi(structFields, candidateNames(i)), 1);
        if ~isempty(idx)
            value = {S.(structFields(idx))}';
            return;
        end
    end
end

function tf = HasFieldOrProperty(obj, name)

    name = char(name);

    if isstruct(obj)
        tf = isfield(obj, name);
    elseif isobject(obj)
        tf = isprop(obj, name);
    else
        tf = false;
    end
end

function value = GetFieldOrProperty(obj, name)

    name = char(name);
    value = obj.(name);
end
