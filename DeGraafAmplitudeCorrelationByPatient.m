function outputs = DeGraafAmplitudeCorrelationByPatient(cfg, varargin)
% DeGraafAmplitudeCorrelationByPatient
%
% Build De Graaf-style amplitude covariance and correlation matrices from
% fitted LCModel spectral basis functions.
%
% For each patient:
%   - find Division_1 part .coord files
%   - load the patient's Division_36_part_1 fitted baseline
%   - append that baseline to each part's fitted spectral basis
%   - invert the complete basis-overlap information matrix
%   - remove the baseline parameter after inversion
%   - parse matching .coord concentration and %SD / CRLB table
%   - mask invalid metabolites
%   - average part-level amplitude correlation matrices
%
% Main outputs:
%
%   outputs.patientResults(1).meanAmplitudeCorrTable
%   outputs.patientResults(1).meanAmplitudeCovTable
%
%   outputs.patientResultsByID.P11.meanAmplitudeCorrTable
%   outputs.patientResultsByID.P11.meanAmplitudeCovTable
%
%   outputs.group.meanAmplitudeCorrTable
%   outputs.group.meanAmplitudeCovTable

    if nargin < 1 || isempty(cfg)
        cfg = ProjectConfig();

    elseif ischar(cfg) || isstring(cfg)
        varargin = [{cfg}, varargin];
        cfg = ProjectConfig();
    end

    opts = ParseDeGraafOptionsFromConfig(cfg, varargin{:});
    patients = BuildPatientFolderTableForDeGraaf(opts);

    metabList = string(opts.metabList(:));
    nMetabs = numel(metabList);
    nPatients = height(patients);

    patientResults = struct([]);
    patientResultsByID = struct;

    patientMeanCorrStack = nan(nMetabs, nMetabs, nPatients);
    patientMeanAbsCorrStack = nan(nMetabs, nMetabs, nPatients);
    patientMeanCovStack = nan(nMetabs, nMetabs, nPatients);

    fprintf('\nDe Graaf / LCModel amplitude correlation analysis\n');
    fprintf('Load mode: %s\n', opts.loadMode);
    fprintf('Number of patient folders: %d\n', nPatients);
    fprintf('Division used: Division_%d\n', opts.division);
    fprintf('Include Division-36 baseline: %d\n', opts.includeDiv36Baseline);
    fprintf('Mask invalid CRLB rows/columns: %d\n', opts.maskInvalidCRLB);

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

        baselineAxis = [];
        baselineVector = [];

        if opts.includeDiv36Baseline
            baselineCoordFile = FindDivisionCoordFile( ...
                coordDir, filePrefix, 36, 1);
            [~, baselineFitData] = VDIIO.ReadLCMCoord(baselineCoordFile);
            [baselineAxis, baselineVector] = ExtractDivisionBaseline( ...
                baselineFitData, baselineCoordFile);
        end

        printFiles = FindDivisionPrintFiles( ...
            coordDir, ...
            filePrefix, ...
            opts.division);

        nParts = numel(printFiles);

        if nParts == 0
            warning('No Division_%d .print files found for patient %s.', ...
                opts.division, patientID);
            continue;
        end

        fprintf('Found %d Division_%d print files.\n', nParts, opts.division);

        partCorrStack = nan(nMetabs, nMetabs, nParts);
        partCovStack = nan(nMetabs, nMetabs, nParts);
        partCRLB = nan(nParts, nMetabs);
        partConc = nan(nParts, nMetabs);
        partSD = nan(nParts, nMetabs);

        partInfo = table( ...
            nan(nParts, 1), ...
            strings(nParts, 1), ...
            strings(nParts, 1), ...
            false(nParts, 1), ...
            strings(nParts, 1), ...
            'VariableNames', { ...
            'part', ...
            'printFile', ...
            'coordFile', ...
            'wasParsed', ...
            'errorMessage'});

    
        skippedParts = table(strings(0,1), strings(0,1), ...
            'VariableNames', {'printFile', 'reason'});

        for partIdx = 1:nParts

            printFile = printFiles(partIdx);
            coordFile = replace(printFile, ".print", ".coord");

            partNumber = ParsePartNumberFromFilename(printFile);

            try
                coordCRLBTable = ParseLCModelCoordCRLB(coordFile);

                [corrMatrix, covMatrix, concVec, crlbVec, sdVec] = ...
                    BuildBasisOverlapAmplitudeMatrices( ...
                        coordFile, ...
                        coordCRLBTable, ...
                        metabList, ...
                        baselineAxis, ...
                        baselineVector, ...
                        opts.includeDiv36Baseline, ...
                        opts.maskInvalidCRLB, ...
                        opts.invalidCRLBValue);

                partCorrStack(:, :, partIdx) = corrMatrix;
                partCovStack(:, :, partIdx) = covMatrix;
                partConc(partIdx, :) = concVec;
                partCRLB(partIdx, :) = crlbVec;
                partSD(partIdx, :) = sdVec;

                partInfo.part(partIdx) = partNumber;
                partInfo.printFile(partIdx) = printFile;
                partInfo.coordFile(partIdx) = coordFile;
                partInfo.wasParsed(partIdx) = true;
                partInfo.errorMessage(partIdx) = "";

            catch ME

                warning('Skipping part file: %s\nReason: %s', printFile, ME.message);

                partInfo.part(partIdx) = partNumber;
                partInfo.printFile(partIdx) = printFile;
                partInfo.coordFile(partIdx) = replace(printFile, ".print", ".coord");
                partInfo.wasParsed(partIdx) = false;
                partInfo.errorMessage(partIdx) = string(ME.message);

                skippedParts = [skippedParts; ...
                    table(printFile, string(ME.message), ...
                    'VariableNames', {'printFile', 'reason'})]; %#ok<AGROW>
            end
        end

        meanAmplitudeCorrMatrix = AverageCorrelationMatricesFisher(partCorrStack);

        % This is mean(abs(part amplitude correlations)), not abs(mean correlation).
        % It captures correlation magnitude even when signs differ between parts.
        meanAbsAmplitudeCorrMatrix = mean(abs(partCorrStack), 3, 'omitnan');

        meanAmplitudeCovMatrix = mean(partCovStack, 3, 'omitnan');

        nValidPartCorrMatrix = sum(~isnan(partCorrStack), 3);
        nValidPartAbsCorrMatrix = sum(~isnan(partCorrStack), 3);
        nValidPartCovMatrix = sum(~isnan(partCovStack), 3);

        crlbSummaryTable = SummarizePartCRLBs( ...
            metabList, ...
            partConc, ...
            partCRLB, ...
            partSD, ...
            opts.invalidCRLBValue);

        crlbLongTable = BuildCRLBLongTableFromPartArrays( ...
            patientID, ...
            opts.division, ...
            metabList, ...
            partInfo, ...
            partConc, ...
            partCRLB, ...
            partSD);

        meanAmplitudeCorrTable = MatrixToMetabTable(meanAmplitudeCorrMatrix, metabList);
        meanAbsAmplitudeCorrTable = MatrixToMetabTable(meanAbsAmplitudeCorrMatrix, metabList);
        meanAmplitudeCovTable = MatrixToMetabTable(meanAmplitudeCovMatrix, metabList);
        nValidPartCorrTable = MatrixToMetabTable(nValidPartCorrMatrix, metabList);
        nValidPartAbsCorrTable = MatrixToMetabTable(nValidPartAbsCorrMatrix, metabList);
        nValidPartCovTable = MatrixToMetabTable(nValidPartCovMatrix, metabList);

        patientResults(pIdx).patientID = patientID;
        patientResults(pIdx).coordDir = coordDir;
        patientResults(pIdx).division = opts.division;
        patientResults(pIdx).metabList = metabList;

        patientResults(pIdx).printFiles = printFiles;
        patientResults(pIdx).partInfo = partInfo;
        patientResults(pIdx).skippedParts = skippedParts;

        patientResults(pIdx).partCorrStack = partCorrStack;
        patientResults(pIdx).partCovStack = partCovStack;
        patientResults(pIdx).partConc = partConc;
        patientResults(pIdx).partCRLB = partCRLB;
        patientResults(pIdx).partSD = partSD;

        patientResults(pIdx).meanAmplitudeCorrMatrix = meanAmplitudeCorrMatrix;
        patientResults(pIdx).meanAbsAmplitudeCorrMatrix = meanAbsAmplitudeCorrMatrix;
        patientResults(pIdx).meanAmplitudeCovMatrix = meanAmplitudeCovMatrix;

        patientResults(pIdx).nValidPartCorrMatrix = nValidPartCorrMatrix;
        patientResults(pIdx).nValidPartAbsCorrMatrix = nValidPartAbsCorrMatrix;
        patientResults(pIdx).nValidPartCovMatrix = nValidPartCovMatrix;

        patientResults(pIdx).meanAmplitudeCorrTable = meanAmplitudeCorrTable;
        patientResults(pIdx).meanAbsAmplitudeCorrTable = meanAbsAmplitudeCorrTable;
        patientResults(pIdx).meanAmplitudeCovTable = meanAmplitudeCovTable;
        patientResults(pIdx).nValidPartCorrTable = nValidPartCorrTable;
        patientResults(pIdx).nValidPartAbsCorrTable = nValidPartAbsCorrTable;
        patientResults(pIdx).nValidPartCovTable = nValidPartCovTable;
        patientResults(pIdx).crlbSummaryTable = crlbSummaryTable;
        patientResults(pIdx).crlbLongTable = crlbLongTable;

        safeField = matlab.lang.makeValidName(char(patientID));
        patientResultsByID.(safeField) = patientResults(pIdx);

        patientMeanCorrStack(:, :, pIdx) = meanAmplitudeCorrMatrix;
        patientMeanAbsCorrStack(:, :, pIdx) = meanAbsAmplitudeCorrMatrix;
        patientMeanCovStack(:, :, pIdx) = meanAmplitudeCovMatrix;
    end

    groupMeanAmplitudeCorrMatrix = AverageCorrelationMatricesFisher(patientMeanCorrStack);

    % This is mean(abs(patient-level amplitude correlations)), not abs(mean correlation).
    groupMeanAbsAmplitudeCorrMatrix = mean(patientMeanAbsCorrStack, 3, 'omitnan');

    groupMeanAmplitudeCovMatrix = mean(patientMeanCovStack, 3, 'omitnan');

    nPatientsCorrMatrix = sum(~isnan(patientMeanCorrStack), 3);
    nPatientsAbsCorrMatrix = sum(~isnan(patientMeanAbsCorrStack), 3);
    nPatientsCovMatrix = sum(~isnan(patientMeanCovStack), 3);

    outputs = struct;
    outputs.settings = opts;
    outputs.patients = patients;
    outputs.metabList = metabList;

    outputs.patientResults = patientResults;
    outputs.patientResultsByID = patientResultsByID;

    outputs.group.patientMeanCorrStack = patientMeanCorrStack;
    outputs.group.patientMeanAbsCorrStack = patientMeanAbsCorrStack;
    outputs.group.patientMeanCovStack = patientMeanCovStack;

    outputs.group.meanAmplitudeCorrMatrix = groupMeanAmplitudeCorrMatrix;
    outputs.group.meanAbsAmplitudeCorrMatrix = groupMeanAbsAmplitudeCorrMatrix;
    outputs.group.meanAmplitudeCovMatrix = groupMeanAmplitudeCovMatrix;

    outputs.group.nPatientsCorrMatrix = nPatientsCorrMatrix;
    outputs.group.nPatientsAbsCorrMatrix = nPatientsAbsCorrMatrix;
    outputs.group.nPatientsCovMatrix = nPatientsCovMatrix;

    outputs.group.meanAmplitudeCorrTable = MatrixToMetabTable(groupMeanAmplitudeCorrMatrix, metabList);
    outputs.group.meanAbsAmplitudeCorrTable = MatrixToMetabTable(groupMeanAbsAmplitudeCorrMatrix, metabList);
    outputs.group.meanAmplitudeCovTable = MatrixToMetabTable(groupMeanAmplitudeCovMatrix, metabList);

    outputs.group.nPatientsCorrTable = MatrixToMetabTable(nPatientsCorrMatrix, metabList);
    outputs.group.nPatientsAbsCorrTable = MatrixToMetabTable(nPatientsAbsCorrMatrix, metabList);
    outputs.group.nPatientsCovTable = MatrixToMetabTable(nPatientsCovMatrix, metabList);
end


function opts = ParseDeGraafOptionsFromConfig(cfg, varargin)

    opts = struct;

    opts.coordRoot = "";
    opts.loadMode = "allSubfolders";

    opts.selectedPatientID = "";
    opts.selectedPatientIDs = strings(0, 1);

    opts.selectedSubfolder = "";
    opts.selectedSubfolders = strings(0, 1);

    opts.division = 1;
    opts.includeDiv36Baseline = true;

    opts.addPatientPrefixToFilenames = true;

    opts.maskInvalidCRLB = true;
    opts.invalidCRLBValue = 999;

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

    % Use cfg.load by default, so this script behaves like Load_subsets_multi_patient_input.
    if HasFieldOrProperty(cfg, "load")
        cfgLoad = GetFieldOrProperty(cfg, "load");

        if HasFieldOrProperty(cfgLoad, "coordDir")
            opts.coordRoot = string(GetFieldOrProperty(cfgLoad, "coordDir"));
        end

        if HasFieldOrProperty(cfgLoad, "mode")
            opts.loadMode = string(GetFieldOrProperty(cfgLoad, "mode"));
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

    % Optional dedicated section:
    % cfg.degraaf.loadMode = ...
    % cfg.degraaf.selectedPatientIDs = ...
    % cfg.degraaf.division = 1
    if HasFieldOrProperty(cfg, "degraaf")
        cfgD = GetFieldOrProperty(cfg, "degraaf");

        if HasFieldOrProperty(cfgD, "coordRoot")
            opts.coordRoot = string(GetFieldOrProperty(cfgD, "coordRoot"));
        end

        if HasFieldOrProperty(cfgD, "coordDir")
            opts.coordRoot = string(GetFieldOrProperty(cfgD, "coordDir"));
        end

        if HasFieldOrProperty(cfgD, "loadMode")
            opts.loadMode = string(GetFieldOrProperty(cfgD, "loadMode"));
        end

        if HasFieldOrProperty(cfgD, "selectedPatientID")
            opts.selectedPatientID = string(GetFieldOrProperty(cfgD, "selectedPatientID"));
        end

        if HasFieldOrProperty(cfgD, "selectedPatientIDs")
            opts.selectedPatientIDs = string(GetFieldOrProperty(cfgD, "selectedPatientIDs"));
        end

        if HasFieldOrProperty(cfgD, "selectedSubfolder")
            opts.selectedSubfolder = string(GetFieldOrProperty(cfgD, "selectedSubfolder"));
        end

        if HasFieldOrProperty(cfgD, "selectedSubfolders")
            opts.selectedSubfolders = string(GetFieldOrProperty(cfgD, "selectedSubfolders"));
        end

        if HasFieldOrProperty(cfgD, "division")
            opts.division = double(GetFieldOrProperty(cfgD, "division"));
        end

        if HasFieldOrProperty(cfgD, "includeDiv36Baseline")
            opts.includeDiv36Baseline = logical(GetFieldOrProperty( ...
                cfgD, "includeDiv36Baseline"));
        end

        if HasFieldOrProperty(cfgD, "metabList")
            opts.metabList = string(GetFieldOrProperty(cfgD, "metabList"));
        end

        if HasFieldOrProperty(cfgD, "maskInvalidCRLB")
            opts.maskInvalidCRLB = logical(GetFieldOrProperty(cfgD, "maskInvalidCRLB"));
        end

        if HasFieldOrProperty(cfgD, "invalidCRLBValue")
            opts.invalidCRLBValue = double(GetFieldOrProperty(cfgD, "invalidCRLBValue"));
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

    p.addParameter('includeDiv36Baseline', opts.includeDiv36Baseline, ...
        @(x) islogical(x) && isscalar(x));

    p.addParameter('addPatientPrefixToFilenames', opts.addPatientPrefixToFilenames, ...
        @(x) islogical(x) && isscalar(x));

    p.addParameter('maskInvalidCRLB', opts.maskInvalidCRLB, ...
        @(x) islogical(x) && isscalar(x));

    p.addParameter('invalidCRLBValue', opts.invalidCRLBValue, ...
        @(x) isnumeric(x) && isscalar(x));

    p.addParameter('metabList', opts.metabList, ...
        @(x) ischar(x) || isstring(x) || iscellstr(x));

    parse(p, varargin{:});

    opts.coordRoot = string(p.Results.coordRoot);
    opts.loadMode = string(p.Results.loadMode);

    opts.selectedPatientID = string(p.Results.selectedPatientID);
    opts.selectedPatientIDs = string(p.Results.selectedPatientIDs);

    opts.selectedSubfolder = string(p.Results.selectedSubfolder);
    opts.selectedSubfolders = string(p.Results.selectedSubfolders);

    opts.division = double(p.Results.division);
    opts.includeDiv36Baseline = logical(p.Results.includeDiv36Baseline);
    opts.addPatientPrefixToFilenames = logical(p.Results.addPatientPrefixToFilenames);

    opts.maskInvalidCRLB = logical(p.Results.maskInvalidCRLB);
    opts.invalidCRLBValue = double(p.Results.invalidCRLBValue);

    opts.metabList = string(p.Results.metabList);
    opts.metabList = opts.metabList(:);

    if strlength(opts.coordRoot) == 0
        error('coordRoot is empty. Set cfg.paths.coordDir, cfg.load.coordDir, or cfg.degraaf.coordRoot.');
    end
end


function patients = BuildPatientFolderTableForDeGraaf(opts)

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
                printInfo = dir(fullfile(curFolder, "*.print"));

                if isempty(printInfo)
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

            printInfo = dir(fullfile(curFolder, "*.print"));

            if isempty(printInfo)
                error('No .print files found in selected folder: %s', curFolder);
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

                printInfo = dir(fullfile(curFolder, "*.print"));

                if isempty(printInfo)
                    warning('Skipping folder with no .print files: %s', curFolder);
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

                printInfo = dir(fullfile(curFolder, "*.print"));

                if isempty(printInfo)
                    warning('Skipping folder with no .print files: %s', curFolder);
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


function printFiles = FindDivisionPrintFiles(coordDir, filePrefix, divisionNumber)

    coordDir = string(coordDir);
    filePrefix = string(filePrefix);

    pattern = filePrefix + "Division_" + string(divisionNumber) + "_part_*.basis.print";
    info = dir(fullfile(coordDir, pattern));

    if isempty(info)
        pattern = filePrefix + "Division_" + string(divisionNumber) + "_*.basis.print";
        info = dir(fullfile(coordDir, pattern));
    end

    if isempty(info)
        printFiles = strings(0, 1);
        return;
    end

    fileNames = string({info.name});
    partNumbers = nan(numel(fileNames), 1);

    for k = 1:numel(fileNames)
        partNumbers(k) = ParsePartNumberFromFilename(fileNames(k));
    end

    [~, order] = sort(partNumbers);
    info = info(order);

    printFiles = strings(numel(info), 1);

    for k = 1:numel(info)
        printFiles(k) = string(fullfile(info(k).folder, info(k).name));
    end
end


function coordFile = FindDivisionCoordFile(coordDir, filePrefix, divisionNumber, partNumber)

    coordDir = string(coordDir);
    filePrefix = string(filePrefix);

    candidateNames = [ ...
        filePrefix + "Division_" + string(divisionNumber) + ...
            "_part_" + string(partNumber) + ".basis.coord", ...
        filePrefix + "Division_" + string(divisionNumber) + ...
            "_" + string(partNumber) + ".basis.coord"];

    found = strings(0, 1);

    for k = 1:numel(candidateNames)
        candidate = string(fullfile(coordDir, candidateNames(k)));
        if isfile(candidate)
            found(end+1, 1) = candidate; %#ok<AGROW>
        end
    end

    found = unique(found, 'stable');

    if numel(found) ~= 1
        error(['Expected exactly one Division_%d part %d .coord file in %s ' ...
            'with prefix "%s"; found %d.'], ...
            divisionNumber, partNumber, coordDir, filePrefix, numel(found));
    end

    coordFile = found(1);
end


function [axisPPM, baseline] = ExtractDivisionBaseline(fitData, coordFile)

    if ~isstruct(fitData) || ~isscalar(fitData) || ...
            ~isfield(fitData, 'axis') || ~isfield(fitData, 'baseline')
        error('Could not read an LCModel axis and baseline from: %s', coordFile);
    end

    axisPPM = double(fitData.axis(:));
    baseline = double(fitData.baseline(:));

    if isempty(axisPPM) || numel(axisPPM) ~= numel(baseline)
        error('Division-36 baseline and ppm axis sizes do not match in: %s', coordFile);
    end

    if ~isreal(baseline) || any(~isfinite(baseline)) || ...
            any(~isfinite(axisPPM)) || norm(baseline) == 0
        error('Division-36 baseline is not a finite, nonzero real spectrum in: %s', ...
            coordFile);
    end

    axisSteps = diff(axisPPM);
    if ~(all(axisSteps > 0) || all(axisSteps < 0))
        error('Division-36 ppm axis is not strictly ordered in: %s', coordFile);
    end
end


function partNumber = ParsePartNumberFromFilename(fileName)

    fileName = string(fileName);

    tok = regexp(fileName, ...
        'Division_\d+_(?:part_)?(\d+)\.basis\.(?:print|coord)$', ...
        'tokens', 'once');

    if isempty(tok)
        partNumber = NaN;
    else
        partNumber = str2double(tok{1});
    end
end


function [corrMatrix, names] = ParseLCModelPrintCorrelation(printFile)

    printFile = string(printFile);

    if ~isfile(printFile)
        error('Print file does not exist: %s', printFile);
    end

    txt = fileread(printFile);
    lines = splitlines(string(txt));

    startIdx = find(contains(lines, "Correlation coefficients"), 1, 'first');

    if isempty(startIdx)
        error('Could not find "Correlation coefficients" section in: %s', printFile);
    end

    sectionLines = lines(startIdx+1:end);

    % First header token gives the first parameter name, usually NAA.
    firstHeaderToken = "";

    for k = 1:numel(sectionLines)
        curLine = strtrim(sectionLines(k));

        if strlength(curLine) == 0
            continue;
        end

        if contains(curLine, "1000Shift") || contains(curLine, "CONC")
            break;
        end

        % Header lines have names but no decimal correlation numbers.
        decimalNums = regexp(curLine, '[-+]?\d+\.\d+(?:[Ee][-+]?\d+)?', 'match');

        if isempty(decimalNums)
            tokens = split(curLine);
            tokens = tokens(strlength(tokens) > 0);

            if ~isempty(tokens)
                firstHeaderToken = tokens(1);
                break;
            end
        end
    end

    if strlength(firstHeaderToken) == 0
        firstHeaderToken = "NAA";
    end

    rowNames = strings(0, 1);
    rowValues = {};

    curName = "";
    curVals = [];

    numberPattern = '[-+]?\d+\.\d+(?:[Ee][-+]?\d+)?';

    for k = 1:numel(sectionLines)

        line = sectionLines(k);

        if contains(line, "1000Shift") || contains(line, "CONC")
            break;
        end

        nums = regexp(line, numberPattern, 'match');
        starts = regexp(line, numberPattern, 'start');

        if isempty(nums)
            continue;
        end

        vals = str2double(nums);

        beforeFirstNumber = extractBefore(line, starts(1));
        label = strtrim(beforeFirstNumber);

        if strlength(label) > 0

            if strlength(curName) > 0
                rowNames(end+1, 1) = curName; %#ok<AGROW>
                rowValues{end+1, 1} = curVals; %#ok<AGROW>
            end

            curName = label;
            curVals = vals(:).';

        else
            curVals = [curVals, vals(:).']; %#ok<AGROW>
        end
    end

    if strlength(curName) > 0
        rowNames(end+1, 1) = curName;
        rowValues{end+1, 1} = curVals;
    end

    names = [firstHeaderToken; rowNames];

    n = numel(names);
    corrMatrix = nan(n, n);

    for i = 1:n
        corrMatrix(i, i) = 1;
    end

    for rowIdx = 2:n

        vals = rowValues{rowIdx - 1};
        expectedCount = rowIdx - 1;

        if numel(vals) < expectedCount
            warning('Correlation row "%s" has fewer values than expected in %s.', ...
                names(rowIdx), printFile);
            expectedCount = numel(vals);
        end

        vals = vals(1:expectedCount);

        corrMatrix(rowIdx, 1:expectedCount) = vals;
        corrMatrix(1:expectedCount, rowIdx) = vals(:);
    end
end


function crlbTable = ParseLCModelCoordCRLB(coordFile)

    coordFile = string(coordFile);

    if ~isfile(coordFile)
        error('Coord file does not exist: %s', coordFile);
    end

    txt = fileread(coordFile);
    lines = splitlines(string(txt));

    headerIdx = find(contains(lines, "Conc.") & ...
                     contains(lines, "%SD") & ...
                     contains(lines, "Metabolite"), 1, 'first');

    if isempty(headerIdx)
        error('Could not find concentration / %%SD table in: %s', coordFile);
    end

    names = strings(0, 1);
    conc = [];
    crlbPercent = [];
    ratioCr = [];

    numberPattern = '[-+]?(?:\d+\.\d*|\.\d+|\d+)(?:[Ee][-+]?\d+)?';

    for k = headerIdx+1:numel(lines)

        line = strtrim(lines(k));

        if strlength(line) == 0
            continue;
        end

        if contains(line, "lines in following misc") || contains(line, "points on ppm-axis")
            break;
        end

        expr = "^(" + numberPattern + ")\s+(\d+)%\s+(" + numberPattern + ")\s+(.+)$";
        tok = regexp(line, expr, 'tokens', 'once');

        if isempty(tok)
            continue;
        end

        names(end+1, 1) = strtrim(string(tok{4})); %#ok<AGROW>
        conc(end+1, 1) = str2double(tok{1}); %#ok<AGROW>
        crlbPercent(end+1, 1) = str2double(tok{2}); %#ok<AGROW>
        ratioCr(end+1, 1) = str2double(tok{3}); %#ok<AGROW>
    end

    if isempty(names)
        error('No metabolite rows parsed from coord file: %s', coordFile);
    end

    sdAbsolute = conc .* crlbPercent ./ 100;

    crlbTable = table(names, conc, crlbPercent, ratioCr, sdAbsolute, ...
        'VariableNames', {'name', 'concentration', 'crlbPercent', 'ratioCr', 'sdAbsolute'});
end


function [corrMatrix, covMatrix, concVec, crlbVec, sdVec] = ...
    BuildBasisOverlapAmplitudeMatrices(coordFile, coordCRLBTable, metabList, ...
    baselineAxis, baselineVector, includeDiv36Baseline, ...
    maskInvalidCRLB, invalidCRLBValue)

    [~, fitData] = VDIIO.ReadLCMCoord(coordFile);

    if ~isstruct(fitData) || ~isscalar(fitData) || ...
            ~isfield(fitData, 'axis') || ~isfield(fitData, 'basisData') || ...
            ~isfield(fitData, 'basisMetName')
        error('Could not read fitted LCModel basis spectra from: %s', coordFile);
    end

    partAxis = double(fitData.axis(:));
    basisData = double(fitData.basisData);
    basisNames = string(fitData.basisMetName(:));

    if size(basisData, 1) ~= numel(partAxis) || ...
            size(basisData, 2) ~= numel(basisNames)
        error('Fitted basis dimensions do not match their ppm axis or names in: %s', ...
            coordFile);
    end

    if ~isreal(basisData) || isempty(basisData) || ...
            any(~isfinite(basisData), 'all')
        error('Fitted basis is not a finite real spectral matrix: %s', coordFile);
    end

    if numel(unique(basisNames)) ~= numel(basisNames)
        error('Fitted basis contains duplicate metabolite names: %s', coordFile);
    end

    covarianceBasis = basisData;
    covarianceNames = basisNames;

    if includeDiv36Baseline
        baselineAxis = double(baselineAxis(:));
        baselineVector = double(baselineVector(:));

        if numel(partAxis) ~= numel(baselineAxis) || ...
                ~isequal(partAxis, baselineAxis)
            error(['Division-36 baseline grid does not exactly match the fitted ' ...
                'metabolite basis grid in: %s'], coordFile);
        end

        if ~isreal(baselineVector) || any(~isfinite(baselineVector))
            error('Division-36 baseline is not a finite real spectrum: %s', ...
                coordFile);
        end

        % EstimateCovarianceMatrix uses raw fitted .coord spectra without
        % column normalization. Keep that convention for the baseline.
        covarianceBasis = [basisData, baselineVector];
        covarianceNames = [basisNames; "Division36Baseline"];
    end

    [~, fullCovariance] = VDILCM.EstimateCovarianceMatrix( ...
        'basisMatrix', covarianceBasis, ...
        'metabNames', covarianceNames);

    nBasis = numel(basisNames);
    expectedCovarianceSize = nBasis + double(includeDiv36Baseline);

    if ~isequal(size(fullCovariance), repmat(expectedCovarianceSize, 1, 2)) || ...
            any(~isfinite(fullCovariance), 'all')
        error('Basis-overlap covariance is invalid for: %s', coordFile);
    end

    % When present, remove the Division-36 nuisance parameter only after the
    % full inverse. Without it, this selects the unchanged full covariance.
    componentScaleCov = fullCovariance(1:nBasis, 1:nBasis);

    metabList = string(metabList(:));
    analysisComponentNames = strings(0, 1);
    for m = 1:numel(metabList)
        analysisComponentNames = [analysisComponentNames; ...
            split(metabList(m), "+")]; %#ok<AGROW>
    end
    analysisComponentNames = unique(analysisComponentNames);
    componentConc = ones(nBasis, 1);
    coordNames = string(coordCRLBTable.name);

    for b = 1:nBasis
        if ismember(basisNames(b), analysisComponentNames)
            idxCoord = find(strcmp(coordNames, basisNames(b)), 1, 'first');
            if isempty(idxCoord) || ...
                    ~isfinite(coordCRLBTable.concentration(idxCoord)) || ...
                    coordCRLBTable.concentration(idxCoord) <= 0
                error(['Missing positive fitted concentration for analysis ' ...
                    'basis component %s in: %s'], basisNames(b), coordFile);
            end
            componentConc(b) = coordCRLBTable.concentration(idxCoord);
        end
    end

    % Exported component spectra include their fitted amplitudes. Transform
    % scale-factor covariance back to concentration-amplitude covariance.
    componentCov = componentScaleCov .* ...
        (componentConc * componentConc.');

    nMetabs = numel(metabList);
    analysisTransform = zeros(nMetabs, nBasis);
    hasBasisRepresentation = false(nMetabs, 1);

    for m = 1:nMetabs
        components = split(metabList(m), "+");
        for c = 1:numel(components)
            idxBasis = find(strcmp(basisNames, components(c)), 1, 'first');
            if ~isempty(idxBasis)
                analysisTransform(m, idxBasis) = 1;
                hasBasisRepresentation(m) = true;
            end
        end
    end

    covMatrix = analysisTransform * componentCov * analysisTransform.';
    covMatrix = (covMatrix + covMatrix.') ./ 2;
    corrMatrix = CovarianceToCorrelation(covMatrix);

    concVec = nan(1, nMetabs);
    crlbVec = nan(1, nMetabs);
    sdVec = nan(1, nMetabs);

    for m = 1:nMetabs
        idxCoord = find(strcmp(coordNames, metabList(m)), 1, 'first');
        if ~isempty(idxCoord)
            concVec(m) = coordCRLBTable.concentration(idxCoord);
            crlbVec(m) = coordCRLBTable.crlbPercent(idxCoord);
            sdVec(m) = coordCRLBTable.sdAbsolute(idxCoord);
        end
    end

    invalid = ~hasBasisRepresentation.' | ...
        isnan(concVec) | isnan(crlbVec) | isnan(sdVec);

    if maskInvalidCRLB
        invalid = invalid | concVec <= 0 | crlbVec >= invalidCRLBValue;
    end

    invalid = invalid | ~isfinite(diag(covMatrix)).' | diag(covMatrix).' <= 0;

    corrMatrix(invalid, :) = NaN;
    corrMatrix(:, invalid) = NaN;
    covMatrix(invalid, :) = NaN;
    covMatrix(:, invalid) = NaN;

    validDiagonal = find(~invalid);
    corrMatrix(sub2ind(size(corrMatrix), validDiagonal, validDiagonal)) = 1;
end


function corrMatrix = CovarianceToCorrelation(covMatrix)

    variances = diag(covMatrix);
    denominator = sqrt(variances * variances.');
    corrMatrix = covMatrix ./ denominator;
    corrMatrix = (corrMatrix + corrMatrix.') ./ 2;
end


function meanCorrMatrix = AverageCorrelationMatricesFisher(corrStack)

    [nRows, nCols, ~] = size(corrStack);

    meanCorrMatrix = nan(nRows, nCols);

    for r = 1:nRows

        for c = 1:nCols

            vals = squeeze(corrStack(r, c, :));
            vals = vals(~isnan(vals));

            if isempty(vals)
                continue;
            end

            if r == c
                meanCorrMatrix(r, c) = 1;
                continue;
            end

            vals(vals >= 1) = 0.999999;
            vals(vals <= -1) = -0.999999;

            zVals = atanh(vals);
            meanCorrMatrix(r, c) = tanh(mean(zVals, 'omitnan'));
        end
    end
end


function summaryTable = SummarizePartCRLBs(metabList, partConc, partCRLB, partSD, maxAllowedCRLB)

    metabList = string(metabList(:));
    nMetabs = numel(metabList);

    nValidParts = nan(nMetabs, 1);
    meanConc = nan(nMetabs, 1);
    medianConc = nan(nMetabs, 1);
    meanCRLB = nan(nMetabs, 1);
    medianCRLB = nan(nMetabs, 1);
    stdCRLB = nan(nMetabs, 1);
    meanSD = nan(nMetabs, 1);

    for m = 1:nMetabs

        concVals = partConc(:, m);
        crlbVals = partCRLB(:, m);
        sdVals = partSD(:, m);

        valid = ~isnan(concVals) & ...
                ~isnan(crlbVals) & ...
                concVals > 0 & ...
                crlbVals < maxAllowedCRLB;

        nValidParts(m) = sum(valid);

        meanConc(m) = mean(concVals(valid), 'omitnan');
        medianConc(m) = median(concVals(valid), 'omitnan');

        meanCRLB(m) = mean(crlbVals(valid), 'omitnan');
        medianCRLB(m) = median(crlbVals(valid), 'omitnan');
        stdCRLB(m) = std(crlbVals(valid), 0, 'omitnan');

        meanSD(m) = mean(sdVals(valid), 'omitnan');
    end

    summaryTable = table( ...
        metabList, ...
        nValidParts, ...
        meanConc, ...
        medianConc, ...
        meanCRLB, ...
        medianCRLB, ...
        stdCRLB, ...
        meanSD, ...
        'VariableNames', { ...
        'metabolite', ...
        'nValidParts', ...
        'meanConc', ...
        'medianConc', ...
        'meanCRLB', ...
        'medianCRLB', ...
        'stdCRLB', ...
        'meanSD'});
end


function crlbLongTable = BuildCRLBLongTableFromPartArrays( ...
    patientID, divisionNumber, metabList, partInfo, partConc, partCRLB, partSD)

    patientID = string(patientID);
    metabList = string(metabList(:));

    nParts = size(partCRLB, 1);
    nMetabs = numel(metabList);

    crlbLongTable = table();

    for partIdx = 1:nParts

        if istable(partInfo) && any(strcmp(partInfo.Properties.VariableNames, "part"))
            partNumber = partInfo.part(partIdx);
        else
            partNumber = partIdx;
        end

        for metabIdx = 1:nMetabs

            newRow = table();

            newRow.patientID = patientID;
            newRow.division = divisionNumber;
            newRow.part = partNumber;
            newRow.metabolite = metabList(metabIdx);
            newRow.concentration = partConc(partIdx, metabIdx);
            newRow.CRLB = partCRLB(partIdx, metabIdx);
            newRow.sdAbsolute = partSD(partIdx, metabIdx);

            crlbLongTable = [crlbLongTable; newRow]; %#ok<AGROW>
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
