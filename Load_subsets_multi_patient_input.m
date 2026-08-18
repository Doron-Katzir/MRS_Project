function outputs = Load_subsets_multi_patient_input(cfg, varargin)

if nargin < 1 || isempty(cfg)
    cfg = ProjectConfig();

elseif ischar(cfg) || isstring(cfg)
    varargin = [{cfg}, varargin];
    cfg = ProjectConfig();
end

runOptions = ParseRunOptionsFromConfig(cfg, varargin{:});
cacheOptions = CoordCacheOptions(cfg, varargin{:});
patients = BuildPatientTableForLoading(cfg, runOptions);

loaderTimer = tic;
cacheStats = InitializeCoordCacheStats(cacheOptions);

fprintf('\nInput mode: %s\n', runOptions.inputMode);
fprintf('Number of patients/files to load: %d\n', height(patients));

analysisMetabList = cfg.metabolites.analysis;
sumMetabList = cfg.metabolites.sum;
metrics = ["sig", "ratioCr", "CRLB"];

referenceDivision = 36;
chosenMetab = "NAA";
ignoreZeros = true;

% Must match the output safety options in Splice_data_multi_patient_safe.m.
% With these settings the loader reads from:
%   <base LCMFit folder>/<patientID>
% and, inside that folder, only reads files starting with:
%   <patientID>_
usePatientSubfolders = runOptions.usePatientSubfolders;
addPatientPrefixToFilenames = runOptions.addPatientPrefixToFilenames;

%% Read and organize each patient separately

patientResults = struct([]);

for pIdx = 1:height(patients)

    patientID = patients.patientID(pIdx);
    baseCoordDir = patients.coordDir(pIdx);
    coordDir = ResolvePatientCoordDir(baseCoordDir, patientID, usePatientSubfolders);

    safePatientID = string(matlab.lang.makeValidName(char(patientID)));
    filePrefix = "";
    if addPatientPrefixToFilenames
        filePrefix = safePatientID + "_";
    end

    fprintf('\nReading patient %s\n', patientID);
    fprintf('Base output dir: %s\n', baseCoordDir);
    fprintf('Resolved patient output dir: %s\n', coordDir);
    fprintf('Filename prefix: %s\n', filePrefix);

    coordStageTimer = tic;
    [quant, fitData, coordFiles, badCoordTable, cacheEvent] = ReadPatientCoordCached( ...
        coordDir, ...
        patientID, ...
        cacheOptions, ...
        'filePrefix', filePrefix);
    cacheEvent.coordStageSeconds = toc(coordStageTimer);
    cacheStats = UpdateCoordCacheStats(cacheStats, cacheEvent);

    fprintf('Coord cache %s: %s', patientID, cacheEvent.status);
    if strlength(cacheEvent.reason) > 0
        fprintf(' - %s', cacheEvent.reason);
    end
    fprintf('\n');
    if ~isempty(badCoordTable)
        warning('Some .coord files were skipped for patient %s.', patientID);
        disp(badCoordTable)
    end
    coordTable = quant.metabTable;

    tablesByMetric = BuildTablesByMetric( ...
        coordTable, ...
        analysisMetabList, ...
        metrics);

    % Single-patient model table for the main concentration metric.
    modelTableSig = BuildConcentrationModelTable( ...
        tablesByMetric.sig, ...
        analysisMetabList, ...
        'ignoreZeros', ignoreZeros);

    [heteroResultsSig, skippedHeteroTable] = FitAllMetabolitesSinglePatient( ...
        modelTableSig, ...
        analysisMetabList, ...
        'referenceDivision', referenceDivision);

    patientResults(pIdx).patientID = patientID; %#ok<SAGROW>
    patientResults(pIdx).baseCoordDir = baseCoordDir;
    patientResults(pIdx).coordDir = coordDir;
    patientResults(pIdx).filePrefix = filePrefix;
    patientResults(pIdx).quant = quant;
    patientResults(pIdx).fitData = fitData;
    patientResults(pIdx).coordFiles = coordFiles;
    patientResults(pIdx).badCoordTable = badCoordTable;
    patientResults(pIdx).coordTable = coordTable;
    patientResults(pIdx).tablesByMetric = tablesByMetric;
    patientResults(pIdx).modelTableSig = modelTableSig;
    patientResults(pIdx).heteroResultsSig = heteroResultsSig;
    patientResults(pIdx).skippedHeteroTable = skippedHeteroTable;
end

%% Single-patient sanity plots for the first patient
% These are useful to verify that the per-patient pipeline still behaves
% like your previous single-patient workflow.
% 
% if ~isempty(patientResults)
% 
%     firstPatient = 1;
% 
%     PlotMetabAcrossDivisions( ...
%         patientResults(firstPatient).tablesByMetric.sig, ...
%         chosenMetab, ...
%         'yLabel', "sig", ...
%         'makeFigure', true);
% 
%     PlotMetabGrandAverageBarsMultiMetric( ...
%         patientResults(firstPatient).tablesByMetric, ...
%         chosenMetab, ...
%         'barWidth', 0.35, ...
%         'yPaddingFrac', 0.30);
% 
%     % Optional spectra plot for the first patient.
%     % Change the division number as needed.
%     PlotDivisionFittedSpectraStack( ...
%         patientResults(firstPatient).fitData, ...
%         patientResults(firstPatient).coordFiles, ...
%         6);
% end

%% Multi-patient bias tables
% Key principle:
%   1. Compute alpha_i inside each patient.
%   2. Combine patient-level alpha_i values across patients.
%   3. Test/summarize across patients, not across subdivision parts.

patientBiasTables = struct;
summaryBiasTables = struct;
skippedMultiTables = struct;

for metricIdx = 1:numel(metrics)

    metricName = metrics(metricIdx);
    metricField = matlab.lang.makeValidName(metricName);

    curPatientTables = GetMetricTablesFromPatientResults( ...
        patientResults, ...
        metricName);

    [patientBiasTables.(metricField), skippedMultiTables.(metricField)] = ...
        ComputeMultiPatientDivisionBias( ...
            curPatientTables, ...
            analysisMetabList, ...
            'patientIDs', patients.patientID, ...
            'referenceDivision', referenceDivision, ...
            'ignoreZeros', ignoreZeros, ...
            'metricName', metricName);

    summaryBiasTables.(metricField) = SummarizeMultiPatientDivisionBias( ...
        patientBiasTables.(metricField), ...
        'referenceDivision', referenceDivision);
end

%% Main multi-patient plots and covariance/correlation matrices

% PlotMultiPatientBiasSpaghetti( ...
%     patientBiasTables.sig, ...
%     chosenMetab, ...
%     'excludeReference', true);
% 
% PlotMultiPatientBiasSummary( ...
%     summaryBiasTables.sig, ...
%     chosenMetab, ...
%     'excludeReference', true);

[covPercentBias_NAA, corrPercentBias_NAA, biasMatrix_NAA] = ...
    ComputeBiasCovarianceForMetabolite( ...
        patientBiasTables.sig, ...
        chosenMetab, ...
        'valueColumn', "percentBias", ...
        'excludeReference', true);

disp('Per-patient bias table for sig:')
disp(patientBiasTables.sig)

disp('Across-patient summary table for sig:')
disp(summaryBiasTables.sig)

disp('Percent-bias covariance matrix for chosen metabolite:')
disp(covPercentBias_NAA)

disp('Percent-bias correlation matrix for chosen metabolite:')
disp(corrPercentBias_NAA)

outputs = struct;
outputs.runOptions = runOptions;
outputs.patients = patients;
outputs.patientResults = patientResults;
outputs.patientBiasTables = patientBiasTables;
outputs.summaryBiasTables = summaryBiasTables;
outputs.skippedMultiTables = skippedMultiTables;
outputs.covPercentBias_NAA = covPercentBias_NAA;
outputs.corrPercentBias_NAA = corrPercentBias_NAA;
outputs.biasMatrix_NAA = biasMatrix_NAA;

% Plot-on-demand handles.
% These allow you to generate plots after loading, from RunExperiment.
outputs.plot.summary = @(metricName, metabName, varargin) ...
    PlotSummaryFromOutputs(outputs, metricName, metabName, varargin{:});

outputs.plot.spaghetti = @(metricName, metabName, varargin) ...
    PlotSpaghettiFromOutputs(outputs, metricName, metabName, varargin{:});

PrintCoordCacheSummary(cacheStats, toc(loaderTimer));
end

%% Functions

function tablesByDivision = coordToMetabs(coordTable, wantedMetabolites, ...
    valueColumn)

    wantedMetabolites = string(wantedMetabolites(:));
    valueColumn = string(valueColumn);

    requiredCols = ["filename", "name", valueColumn];
    tableCols = string(coordTable.Properties.VariableNames);

    for c = requiredCols
        if ~ismember(c, tableCols)
            error('coordTable is missing required column "%s".', c);
        end
    end

    % Work on copy
    T = coordTable;
    T.filename = string(T.filename);
    T.name = string(T.name);

    % Keep only wanted metabolites
    T = T(ismember(T.name, wantedMetabolites), :);

    if isempty(T)
        error(['None of the wanted metabolites were found in ' ...
            'coordTable.name.']);
    end

    % Extract division size and part number from filename.
    % Supports both unprefixed and patient-prefixed filenames:
    %   Division_18_1.basis.coord
    %   Division_18_part_1.basis.coord
    %   P01_Division_18_part_1.basis.coord
    nRows = height(T);

    divisionSize = nan(nRows, 1);
    partNumber = nan(nRows, 1);

    for i = 1:nRows

        curFile = T.filename(i);

        tok = regexp(curFile, ...
            'Division_(\d+)_(?:part_)?(\d+)\.basis\.coord', ...
            'tokens', 'once');

        if isempty(tok)
            error('Could not parse filename: %s', curFile);
        end

        divisionSize(i) = str2double(tok{1});
        partNumber(i) = str2double(tok{2});
    end

    T.divisionSize = divisionSize;
    T.partNumber = partNumber;

    % Create one table per division size
    uniqueDivisions = unique(T.divisionSize);
    uniqueDivisions = sort(uniqueDivisions);

    tablesByDivision = struct;

    for d = 1:numel(uniqueDivisions)

        curDivision = uniqueDivisions(d);

        % Rows belonging only to this division size
        Td = T(T.divisionSize == curDivision, :);

        % Unique parts for this division
        parts = unique(Td.partNumber);
        parts = sort(parts);

        % Initialize output table for this division
        outTable = table;
        outTable.name = wantedMetabolites;

        % Add one column per part
        for p = 1:numel(parts)
            colName = sprintf('part_%d', parts(p));
            outTable.(colName) = nan(numel(wantedMetabolites), 1);
        end

        % Fill values
        for m = 1:numel(wantedMetabolites)

            curMetab = wantedMetabolites(m);

            for p = 1:numel(parts)

                curPart = parts(p);
                colName = sprintf('part_%d', curPart);

                idx = Td.name == curMetab & Td.partNumber == curPart;

                if sum(idx) == 1
                    outTable.(colName)(m) = Td.(valueColumn)(idx);

                elseif sum(idx) > 1
                    warning(['Multiple values found for %s,' ...
                        ' Division %d, part %d. Using first value.'], ...
                        curMetab, curDivision, curPart);

                    vals = Td.(valueColumn)(idx);
                    outTable.(colName)(m) = vals(1);

                else
                    % Missing metabolite in this division/part
                    outTable.(colName)(m) = NaN;
                end
            end
        end

        % Store in struct
        fieldName = sprintf('Division_%d', curDivision);
        fieldName = matlab.lang.makeValidName(fieldName);

        tablesByDivision.(fieldName) = outTable;
    end
end

function PlotSummaryFromOutputs(outputs, metricName, metabName, varargin)

    metricName = string(metricName);
    metricField = matlab.lang.makeValidName(metricName);
    
    if ~isfield(outputs.summaryBiasTables, metricField)
        error('Metric "%s" was not found in outputs.summaryBiasTables.', metricName);
    end
    
    PlotMultiPatientBiasSummary( ...
        outputs.summaryBiasTables.(metricField), ...
        metabName, ...
        'metricName', metricName, ...
        varargin{:});
end

function PlotSpaghettiFromOutputs(outputs, metricName, metabName, varargin)

    metricName = string(metricName);
    metricField = matlab.lang.makeValidName(metricName);
    
    if ~isfield(outputs.patientBiasTables, metricField)
        error('Metric "%s" was not found in outputs.patientBiasTables.', metricName);
    end
    
    PlotMultiPatientBiasSpaghetti( ...
        outputs.patientBiasTables.(metricField), ...
        metabName, ...
        'metricName', metricName, ...
        varargin{:});
end

function summary = PlotMetabAcrossDivisions(tablesByDivisionInput, metabName, varargin)
% PlotMetabAcrossDivisions
%
% Plots one metabolite across all divisions.
%
% Supports:
%   1. Single dataset:
%       PlotMetabAcrossDivisions(tablesByDivision, "NAA")
%
%   2. Multiple datasets / patients in the future:
%       PlotMetabAcrossDivisions({tablesByDivisionPatient1, ...
%                                 tablesByDivisionPatient2}, "NAA")
%
% Each tablesByDivision structure should contain fields like:
%   Division_2
%   Division_4
%   Division_6
%   Division_9
%   Division_12
%   Division_18
%   Division_36
%
% Each field should contain a table like:
%   name    part_1    part_2    part_3 ...
%
% Zeros are shown in raw patient traces but ignored in averages.

    % ------------------------------------------------------------
    % Optional inputs
    % ------------------------------------------------------------
    p = inputParser;
    p.addParameter('timePerSet', 1, @(x) isnumeric(x) && isscalar(x));
    p.addParameter('timeUnit', "sets", @(x) ischar(x) || isstring(x));
    p.addParameter('showPatients', true, @(x) islogical(x) && isscalar(x));
    p.addParameter('showGrandMean', true, @(x) islogical(x) && isscalar(x));
    p.addParameter('yLabel', "sig", @(x) ischar(x) || isstring(x));
    p.addParameter('makeFigure', true, @(x) islogical(x) && isscalar(x));
    p.addParameter('barWidth', 0.45, @(x) isnumeric(x) && isscalar(x));
    p.addParameter('zoomY', true, @(x) islogical(x) && isscalar(x));
    p.addParameter('yPaddingFrac', 0.15, @(x) isnumeric(x) && isscalar(x));
    parse(p, varargin{:});

    timePerSet = p.Results.timePerSet;
    timeUnit = string(p.Results.timeUnit);
    showPatients = p.Results.showPatients;
    showGrandMean = p.Results.showGrandMean;
    yLabelText = string(p.Results.yLabel);
    makeFigure = p.Results.makeFigure;

    metabName = string(metabName);

    % ------------------------------------------------------------
    % Convert input to cell array of patients/datasets
    % ------------------------------------------------------------
    if iscell(tablesByDivisionInput)
        patientTables = tablesByDivisionInput(:);

    elseif isstruct(tablesByDivisionInput) && numel(tablesByDivisionInput) > 1
        patientTables = arrayfun(@(s) s, tablesByDivisionInput, ...
            'UniformOutput', false);

    elseif isstruct(tablesByDivisionInput)
        patientTables = {tablesByDivisionInput};

    else
        error('Input must be a tablesByDivision struct, struct array, or cell array of structs.');
    end

    nPatients = numel(patientTables);

    % ------------------------------------------------------------
    % Collect all division field names across patients
    % ------------------------------------------------------------
    divisionFields = strings(0, 1);

    for pIdx = 1:nPatients
        curFields = string(fieldnames(patientTables{pIdx}));

        for fIdx = 1:numel(curFields)
            if ~any(divisionFields == curFields(fIdx))
                divisionFields(end+1, 1) = curFields(fIdx); %#ok<AGROW>
            end
        end
    end

    % Sort fields numerically by division number
    divisionNumbers = nan(numel(divisionFields), 1);

    for i = 1:numel(divisionFields)
        divisionNumbers(i) = ParseDivisionNumber(divisionFields(i));
    end

    [~, order] = sort(divisionNumbers);
    divisionFields = divisionFields(order);
    divisionNumbers = divisionNumbers(order);

    % ------------------------------------------------------------
    % Prepare figure
    % ------------------------------------------------------------
    if makeFigure
        figure;
        hold on;
    end

    summary = struct;
    summary.metabName = metabName;
    summary.nPatients = nPatients;

    % ------------------------------------------------------------
    % Main loop over divisions
    % ------------------------------------------------------------
    for d = 1:numel(divisionFields)

        curField = divisionFields(d);
        nAvg = divisionNumbers(d);

        if isnan(nAvg)
            warning('Could not parse division number from field "%s". Skipping.', curField);
            continue;
        end

        % --------------------------------------------------------
        % Extract this metabolite from each patient for this division
        % --------------------------------------------------------
        patientVectors = cell(nPatients, 1);

        for pIdx = 1:nPatients

            curStruct = patientTables{pIdx};

            if ~isfield(curStruct, curField)
                patientVectors{pIdx} = [];
                continue;
            end

            curTable = curStruct.(curField);

            y = ExtractMetabVectorFromTable(curTable, metabName);

            patientVectors{pIdx} = y;
        end

        % Remove patients that do not have this division/metabolite
        hasData = cellfun(@(x) ~isempty(x), patientVectors);
        patientVectors = patientVectors(hasData);

        if isempty(patientVectors)
            warning('No data found for metabolite "%s" in %s.', metabName, curField);
            continue;
        end

        % --------------------------------------------------------
        % Convert patient vectors into matrix:
        %   rows    = patients
        %   columns = time parts
        % --------------------------------------------------------
        nUsedPatients = numel(patientVectors);
        nParts = max(cellfun(@numel, patientVectors));

        Yraw = nan(nUsedPatients, nParts);

        for pIdx = 1:nUsedPatients
            curY = patientVectors{pIdx};
            Yraw(pIdx, 1:numel(curY)) = curY;
        end

        % --------------------------------------------------------
        % Build time axis
        %
        % If Division_6 means average 6 original sets per part:
        %   part_1 center = 3
        %   part_2 center = 9
        %   part_3 center = 15
        %
        % Formula:
        %   time = ((part index - 0.5) * nAvg) * timePerSet
        % --------------------------------------------------------
        partIdx = 1:nParts;
        timeAxis = ((partIdx - 0.5) * nAvg) * timePerSet;

        % --------------------------------------------------------
        % Compute averages while ignoring zeros
        % --------------------------------------------------------
        YforMean = Yraw;
        YforMean(YforMean == 0) = NaN;

        meanAcrossPatients = mean(YforMean, 1, 'omitnan');
        grandMean = mean(YforMean(:), 'omitnan');

        nNonZeroPerTimePoint = sum(~isnan(YforMean), 1);

        % --------------------------------------------------------
        % Store summary
        % --------------------------------------------------------
        outField = char(curField);

        summary.(outField).division = nAvg;
        summary.(outField).time = timeAxis;
        summary.(outField).rawValues = Yraw;
        summary.(outField).valuesForMean = YforMean;
        summary.(outField).meanAcrossPatients = meanAcrossPatients;
        summary.(outField).grandMean = grandMean;
        summary.(outField).nNonZeroPerTimePoint = nNonZeroPerTimePoint;
        summary.(outField).nPatientsUsed = nUsedPatients;

        % --------------------------------------------------------
        % Plot
        % --------------------------------------------------------
        if makeFigure

            ax = gca;
            colorList = colororder(ax);
            thisColor = colorList(mod(d-1, size(colorList,1)) + 1, :);
        
            % Raw patient traces, including zeros
            if showPatients
                for pIdx = 1:nUsedPatients
                    plot(timeAxis, Yraw(pIdx, :), ':o', ...
                        'Color', thisColor, ...
                        'HandleVisibility', 'off');
                end
            end
        
            % Average time series, zeros ignored
            plot(timeAxis, meanAcrossPatients, '-o', ...
                'Color', thisColor, ...
                'LineWidth', 2, ...
                'DisplayName', sprintf('%s mean', curField));
        
            % Grand mean across patients and time, zeros ignored
            if showGrandMean && ~isnan(grandMean)
                plot([timeAxis(1), timeAxis(end)], [grandMean, grandMean], '--', ...
                    'Color', thisColor, ...
                    'LineWidth', 1.5, ...
                    'DisplayName', sprintf('%s grand mean', curField));
            end
        end
    end

    % ------------------------------------------------------------
    % Final plot formatting
    % ------------------------------------------------------------
    if makeFigure
        xlabel(sprintf('Time (%s)', timeUnit));
        ylabel(yLabelText, 'Interpreter', 'none');
        title(sprintf('%s across divisions', metabName), 'Interpreter', 'none');
        legend('Location', 'northeast', 'Interpreter', 'none');
        grid on;
        hold off;
    end
end

function nAvg = ParseDivisionNumber(fieldName)
    % Extracts the number from field names like:
    %   Division_2
    %   Division_6
    %   Division_18
    
    fieldName = string(fieldName);
    
    tok = regexp(fieldName, 'Division_(\d+)', 'tokens', 'once');
    
    if isempty(tok)
        nAvg = NaN;
    else
        nAvg = str2double(tok{1});
    end
end

function y = ExtractMetabVectorFromTable(T, metabName)
    % Extracts one metabolite row from a Division_* table.
    %
    % Expected table format:
    %   name    part_1    part_2    part_3 ...
    
    y = [];
    
    metabName = string(metabName);
    
    if ~ismember("name", string(T.Properties.VariableNames))
        error('The table does not contain a column named "name".');
    end
    
    names = string(T.name);
    
    rowIdx = find(strcmpi(names, metabName), 1, 'first');
    
    if isempty(rowIdx)
        warning('Metabolite "%s" was not found in this table.', metabName);
        return;
    end
    
    % Extract all columns except the first one, which is the metabolite name
    y = T{rowIdx, 2:end};
    
    % Make sure output is numeric row vector
    y = double(y);
    y = reshape(y, 1, []);
end

function grandTable = PlotMetabGrandAverageBars(tablesByDivisionInput, ...
    metabList, varargin)
    % PlotMetabGrandAverageBars
    %
    % For each chosen metabolite, plots a bar plot:
    %   x-axis    = division size
    %   bar value = grand average across time, ignoring zeros
    %   error bar = standard deviation across time, ignoring zeros
    %
    % Supports:
    %   Single dataset:
    %       PlotMetabGrandAverageBars(tablesByDivision, "NAA")
    %
    %   Multiple metabolites:
    %       PlotMetabGrandAverageBars(tablesByDivision, ["NAA", "Glu", "Cr+PCr"])
    %
    %   Future multiple patients:
    %       PlotMetabGrandAverageBars({tablesPatient1, tablesPatient2}, "NAA")
    
        p = inputParser;
        p.addParameter('divisions', [], @(x) isempty(x) || isnumeric(x));
        p.addParameter('ignoreZeros', true, @(x) islogical(x) && isscalar(x));
        p.addParameter('yLabel', "sig", @(x) ischar(x) || isstring(x));
        p.addParameter('makeFigure', true, @(x) islogical(x) && isscalar(x));
        parse(p, varargin{:});
    
        requestedDivisions = p.Results.divisions;
        ignoreZeros = p.Results.ignoreZeros;
        yLabelText = string(p.Results.yLabel);
        makeFigure = p.Results.makeFigure;
    
        metabList = string(metabList(:));
    
        % ------------------------------------------------------------
        % Convert input to cell array of patients/datasets
        % ------------------------------------------------------------
        if iscell(tablesByDivisionInput)
            patientTables = tablesByDivisionInput(:);
    
        elseif isstruct(tablesByDivisionInput) && numel(tablesByDivisionInput) > 1
            patientTables = num2cell(tablesByDivisionInput);
    
        elseif isstruct(tablesByDivisionInput)
            patientTables = {tablesByDivisionInput};
    
        else
            error('Input must be a tablesByDivision struct, struct array, or cell array of structs.');
        end
    
        nPatients = numel(patientTables);
    
        % ------------------------------------------------------------
        % Collect all division fields
        % ------------------------------------------------------------
        divisionFields = strings(0, 1);
    
        for pIdx = 1:nPatients
            curFields = string(fieldnames(patientTables{pIdx}));
    
            for fIdx = 1:numel(curFields)
                if ~any(divisionFields == curFields(fIdx))
                    divisionFields(end+1, 1) = curFields(fIdx); %#ok<AGROW>
                end
            end
        end
    
        divisionNumbers = nan(numel(divisionFields), 1);
    
        for i = 1:numel(divisionFields)
            divisionNumbers(i) = ParseDivisionNumber(divisionFields(i));
        end
    
        [~, order] = sort(divisionNumbers);
        divisionFields = divisionFields(order);
        divisionNumbers = divisionNumbers(order);
    
        % Optional: keep only requested divisions
        if ~isempty(requestedDivisions)
            keepIdx = ismember(divisionNumbers, requestedDivisions);
    
            if ~any(keepIdx)
                error('None of the requested divisions were found.');
            end
    
            divisionFields = divisionFields(keepIdx);
            divisionNumbers = divisionNumbers(keepIdx);
        end
    
        % ------------------------------------------------------------
        % Compute grand mean and std
        % ------------------------------------------------------------
        rows = {};
    
        for m = 1:numel(metabList)
    
            metabName = metabList(m);
    
            grandMean = nan(numel(divisionFields), 1);
            grandStd  = nan(numel(divisionFields), 1);
            nValid    = nan(numel(divisionFields), 1);
    
            for d = 1:numel(divisionFields)
    
                curField = divisionFields(d);
                curDivision = divisionNumbers(d);
    
                allValues = [];
    
                for pIdx = 1:nPatients
    
                    curStruct = patientTables{pIdx};
    
                    if ~isfield(curStruct, curField)
                        continue;
                    end
    
                    curTable = curStruct.(curField);
    
                    y = ExtractMetabVectorFromTable(curTable, metabName);
    
                    if isempty(y)
                        continue;
                    end
    
                    y = double(y(:));
    
                    if ignoreZeros
                        y(y == 0) = NaN;
                    end
    
                    allValues = [allValues; y]; %#ok<AGROW>
                end
    
                validValues = allValues(~isnan(allValues));
    
                if ~isempty(validValues)
                    grandMean(d) = mean(validValues);
                    grandStd(d)  = std(validValues, 0);
                    nValid(d)    = numel(validValues);
                else
                    grandMean(d) = NaN;
                    grandStd(d)  = NaN;
                    nValid(d)    = 0;
                end
    
                rows(end+1, :) = { ...
                    metabName, ...
                    curField, ...
                    curDivision, ...
                    grandMean(d), ...
                    grandStd(d), ...
                    nValid(d), ...
                    nPatients}; %#ok<AGROW>
            end
    
            % --------------------------------------------------------
            % Plot one bar figure per metabolite
            % --------------------------------------------------------
            if makeFigure
            
                figure;
                hold on;
            
                x = 1:numel(divisionNumbers);
                xLabels = string(divisionNumbers);
            
                bar(x, grandMean);
            
                errorbar(x, grandMean, grandStd, ...
                    'k', ...
                    'LineStyle', 'none', ...
                    'LineWidth', 1.5);
            
                xticks(x);
                xticklabels(xLabels);
            
                xlabel('Division size / sets averaged');
                ylabel(yLabelText, 'Interpreter', 'none');
            
                title(sprintf('%s grand average across time', metabName), ...
                    'Interpreter', 'none');
            
                grid on;
                hold off;
            end
        end
    
        grandTable = cell2table(rows, ...
            'VariableNames', { ...
            'metabolite', ...
            'divisionField', ...
            'division', ...
            'grandMean', ...
            'grandStd', ...
            'nValidValues', ...
            'nPatients'});
    
        grandTable.metabolite = string(grandTable.metabolite);
        grandTable.divisionField = string(grandTable.divisionField);
end


function grandTablesByMetric = PlotMetabGrandAverageBarsMultiMetric( ...
    tablesByMetric, metabList, varargin)
    % PlotMetabGrandAverageBarsMultiMetric
    %
    % For each metabolite, plots bar plots for:
    %   sig
    %   ratioCr
    %   CRLB
    %
    % Input:
    %   tablesByMetric.sig
    %   tablesByMetric.ratioCr
    %   tablesByMetric.CRLB
    %
    % Each field should be a tablesByDivision struct.
    
    p = inputParser;
    p.addParameter('metrics', ["sig", "ratioCr", "CRLB"], ...
        @(x) ischar(x) || isstring(x) || iscellstr(x));
    p.addParameter('divisions', [], @(x) isempty(x) || isnumeric(x));
    p.addParameter('ignoreZeros', true, @(x) islogical(x) && isscalar(x));
    p.addParameter('barWidth', 0.45, @(x) isnumeric(x) && isscalar(x));
    p.addParameter('zoomY', true, @(x) islogical(x) && isscalar(x));
    p.addParameter('yPaddingFrac', 0.15, @(x) isnumeric(x) && isscalar(x));
    parse(p, varargin{:});
    
    metrics = string(p.Results.metrics);
    requestedDivisions = p.Results.divisions;
    ignoreZeros = p.Results.ignoreZeros;
    barWidth = p.Results.barWidth;
    zoomY = p.Results.zoomY;
    yPaddingFrac = p.Results.yPaddingFrac;
    
    metabList = string(metabList(:));
    
    grandTablesByMetric = struct;
    
    for m = 1:numel(metabList)
    
        metabName = metabList(m);
    
        figure;
        tiledlayout(numel(metrics), 1, 'TileSpacing', 'compact');
    
        for k = 1:numel(metrics)
    
            metricName = metrics(k);
            metricField = matlab.lang.makeValidName(metricName);
    
            if ~isfield(tablesByMetric, metricField)
                error('tablesByMetric is missing field "%s".', metricField);
            end
    
            nexttile;
    
            grandTable = PlotMetabGrandAverageBars( ...
                tablesByMetric.(metricField), ...
                metabName, ...
                'divisions', requestedDivisions, ...
                'ignoreZeros', ignoreZeros, ...
                'yLabel', metricName, ...
                'makeFigure', false);
    
            % Store table
            metabField = matlab.lang.makeValidName(metabName);
    
            if ~isfield(grandTablesByMetric, metricField)
                grandTablesByMetric.(metricField) = struct;
            end
    
            grandTablesByMetric.(metricField).(metabField) = grandTable;
    
            % Plot manually into current tile
            Tm = grandTable(grandTable.metabolite == metabName, :);
            Tm = sortrows(Tm, 'division');
    
            x = 1:height(Tm);
            xLabels = string(Tm.division);
    
            bar(x, Tm.grandMean, barWidth);
            hold on;
            
            errorbar(x, Tm.grandMean, Tm.grandStd, ...
                'k', ...
                'LineStyle', 'none', ...
                'LineWidth', 1.5);
            
            xticks(x);
            xticklabels(xLabels);
            
            xlabel('Division size / sets averaged');
            ylabel(metricName, 'Interpreter', 'none');
            
            title(sprintf('%s: %s', metabName, metricName), ...
                'Interpreter', 'none');
            
            grid on;
            
            % ------------------------------------------------------------
            % Zoom y-axis around the data so small differences are visible
            % ------------------------------------------------------------
            if zoomY
            
                yLow = Tm.grandMean - Tm.grandStd;
                yHigh = Tm.grandMean + Tm.grandStd;
            
                yLow = yLow(~isnan(yLow));
                yHigh = yHigh(~isnan(yHigh));
            
                if ~isempty(yLow) && ~isempty(yHigh)
            
                    yMin = min(yLow);
                    yMax = max(yHigh);
            
                    yRange = yMax - yMin;
            
                    if yRange == 0
                        % Handles almost-flat data, e.g. ratioCr for Cr+PCr
                        yRange = abs(yMax) * 0.05;
            
                        if yRange == 0
                            yRange = 1;
                        end
                    end
            
                    pad = yPaddingFrac * yRange;
            
                    ylim([yMin - pad, yMax + pad]);
                end
            end
            
            hold off;
        end
    end
end

function PlotDivisionFittedSpectraStack(fitData, coordFiles, divisionNumber, varargin)
% PlotDivisionFittedSpectraStack
%
% LCModel-style stacked plot where each line is one fitted spectrum
% from a subdivision part.
%
% Example:
%   PlotDivisionFittedSpectraStack(fitData, coordFiles, 6);
%
% Supports patient-prefixed filenames such as:
%   P01_Division_6_part_1.basis.coord
%
% Optional:
%   'dataField'        - "fitData", "phasedData", "residual", "baseline"
%   'ppmRange'         - [0.2 4.2]
%   'normalize'        - true/false
%   'spacingFactor'    - controls vertical spacing
%   'reverseXAxis'     - true/false
%   'lineWidth'        - curve line width

    p = inputParser;
    p.addParameter('dataField', "fitData", @(x) ischar(x) || isstring(x));
    p.addParameter('ppmRange', [0.2 4.2], ...
        @(x) isempty(x) || (isnumeric(x) && numel(x) == 2));
    p.addParameter('normalize', false, @(x) islogical(x) && isscalar(x));
    p.addParameter('spacingFactor', 1.15, @(x) isnumeric(x) && isscalar(x));
    p.addParameter('reverseXAxis', true, @(x) islogical(x) && isscalar(x));
    p.addParameter('lineWidth', 1.0, @(x) isnumeric(x) && isscalar(x));
    parse(p, varargin{:});

    dataField = string(p.Results.dataField);
    ppmRange = p.Results.ppmRange;
    normalize = p.Results.normalize;
    spacingFactor = p.Results.spacingFactor;
    reverseXAxis = p.Results.reverseXAxis;
    lineWidth = p.Results.lineWidth;

    coordFiles = string(coordFiles(:));

    if numel(coordFiles) ~= numel(fitData)
        error('coordFiles and fitData must have the same number of elements.');
    end

    if ~isfield(fitData, char(dataField))
        error('fitData does not contain field "%s". Available fields are: %s', ...
            dataField, strjoin(string(fieldnames(fitData)), ", "));
    end

    % ------------------------------------------------------------
    % Parse division and part from filenames
    % ------------------------------------------------------------
    nFiles = numel(coordFiles);
    division = nan(nFiles, 1);
    part = nan(nFiles, 1);

    for i = 1:nFiles

        [~, baseName, ext] = fileparts(coordFiles(i));
        curName = baseName + ext;

        tok = regexp(curName, ...
            'Division_(\d+)_(?:part_)?(\d+)\.basis\.coord', ...
            'tokens', 'once');

        if isempty(tok)
            warning('Could not parse filename: %s', curName);
            continue;
        end

        division(i) = str2double(tok{1});
        part(i) = str2double(tok{2});
    end

    % ------------------------------------------------------------
    % Select requested division
    % ------------------------------------------------------------
    idx = find(division == divisionNumber);

    if isempty(idx)
        error('No fitted spectra found for Division_%d.', divisionNumber);
    end

    [~, order] = sort(part(idx));
    idx = idx(order);
    selectedParts = part(idx);

    nParts = numel(idx);

    % ------------------------------------------------------------
    % First pass: collect spectra
    % ------------------------------------------------------------
    spectra = cell(nParts, 1);
    xCell = cell(nParts, 1);
    maxAbsVals = nan(nParts, 1);

    for k = 1:nParts

        curIdx = idx(k);

        x = double(fitData(curIdx).axis(:));
        y = real(double(fitData(curIdx).(char(dataField))(:)));

        if ~isempty(ppmRange)
            lo = min(ppmRange);
            hi = max(ppmRange);
            keep = x >= lo & x <= hi;

            x = x(keep);
            y = y(keep);
        end

        if normalize
            scaleVal = max(abs(y), [], 'omitnan');

            if scaleVal ~= 0 && ~isnan(scaleVal)
                y = y ./ scaleVal;
            end
        end

        xCell{k} = x;
        spectra{k} = y;
        maxAbsVals(k) = max(abs(y), [], 'omitnan');
    end

    % Vertical spacing
    baseSpacing = max(maxAbsVals, [], 'omitnan') * spacingFactor;

    if baseSpacing == 0 || isnan(baseSpacing)
        baseSpacing = 1;
    end

    % ------------------------------------------------------------
    % Plot
    % ------------------------------------------------------------
    figure;
    hold on;

    colorMap = lines(nParts);

    for k = 1:nParts

        x = xCell{k};
        y = spectra{k};

        % part_1 at top, later parts below
        yOffset = -(k - 1) * baseSpacing;

        plot(x, y + yOffset, ...
            'Color', colorMap(k, :), ...
            'LineWidth', lineWidth);

        % Label on right side of plot
        labelX = min(x);

        text(labelX, yOffset, sprintf('part_%d', selectedParts(k)), ...
            'HorizontalAlignment', 'left', ...
            'VerticalAlignment', 'middle', ...
            'Interpreter', 'none', ...
            'Color', colorMap(k, :));
    end

    if reverseXAxis
        set(gca, 'XDir', 'reverse');
    end

    xlabel('ppm');

    if normalize
        ylabel(sprintf('%s, normalized and stacked', dataField), ...
            'Interpreter', 'none');
    else
        ylabel(sprintf('%s, stacked', dataField), ...
            'Interpreter', 'none');
    end

    title(sprintf('Division_%d stacked fitted spectra', divisionNumber), ...
        'Interpreter', 'none');

    grid on;
    box on;

    % Remove y tick labels because vertical position is artificial
    yticks([]);

    hold off;
end

function modelTable = BuildConcentrationModelTable(tablesByDivision, metabList, varargin)
% BuildConcentrationModelTable
%
% Converts tablesByDivision into long format.
%
% Each row is one LCModel concentration value, X_ij:
%
%   metabolite | division | divisionID | part | concentration
%
% This table is the input for:
%
%   X_ij = beta0 + alpha_i + epsilon_ij

    p = inputParser;
    p.addParameter('ignoreZeros', true, @(x) islogical(x) && isscalar(x));
    parse(p, varargin{:});

    ignoreZeros = p.Results.ignoreZeros;
    metabList = string(metabList(:));

    divisionFields = string(fieldnames(tablesByDivision));

    rows = {};

    for d = 1:numel(divisionFields)

        curField = divisionFields(d);
        curDivision = ParseDivisionNumber(curField);

        if isnan(curDivision)
            continue;
        end

        curTable = tablesByDivision.(curField);

        for m = 1:numel(metabList)

            metabName = metabList(m);
            values = ExtractMetabVectorFromTable(curTable, metabName);

            if isempty(values)
                continue;
            end

            values = double(values(:));

            for partIdx = 1:numel(values)

                curValue = values(partIdx);

                if ignoreZeros && curValue == 0
                    curValue = NaN;
                end

                if isnan(curValue)
                    continue;
                end

                rows(end+1, :) = { ...
                    metabName, ...
                    curDivision, ...
                    "Division_" + string(curDivision), ...
                    partIdx, ...
                    curValue}; %#ok<AGROW>
            end
        end
    end

    modelTable = cell2table(rows, ...
        'VariableNames', { ...
        'metabolite', ...
        'division', ...
        'divisionID', ...
        'part', ...
        'concentration'});

    modelTable.metabolite = string(modelTable.metabolite);
    modelTable.divisionID = categorical(string(modelTable.divisionID));
end

function result = FitDivisionHeteroscedasticFixedEffect(modelTable, metabName, varargin)
% FitDivisionHeteroscedasticFixedEffect
%
% Implements your mentor's model:
%
%   X_ij = beta0 + alpha_i + epsilon_ij
%
% with group-specific residual variance:
%
%   epsilon_ij ~ N(0, sigma_i^2)
%
% where:
%   i = division
%   j = part within division
%
% beta0 is the reference division mean.
% alpha_i is the bias of division i relative to the reference.
% sigma_i^2 is estimated separately for each division.

    p = inputParser;
    p.addParameter('referenceDivision', 36, @(x) isnumeric(x) && isscalar(x));
    parse(p, varargin{:});

    referenceDivision = p.Results.referenceDivision;
    metabName = string(metabName);

    % Keep only this metabolite.
    Tm = modelTable(strcmpi(modelTable.metabolite, metabName), :);

    if isempty(Tm)
        error('Metabolite "%s" was not found in modelTable.', metabName);
    end

    Tm = sortrows(Tm, {'division', 'part'});

    divisions = unique(Tm.division);
    divisions = sort(divisions(:));

    if ~ismember(referenceDivision, divisions)
        error('Reference division %d was not found for metabolite "%s".', ...
            referenceDivision, metabName);
    end

    nDivisions = numel(divisions);

    nValues = nan(nDivisions, 1);
    meanValue = nan(nDivisions, 1);
    stdValue = nan(nDivisions, 1);
    varValue = nan(nDivisions, 1);
    semValue = nan(nDivisions, 1);

    % ------------------------------------------------------------
    % Estimate one mean and one variance per division.
    % ------------------------------------------------------------
    for i = 1:nDivisions

        curDivision = divisions(i);

        curValues = Tm.concentration(Tm.division == curDivision);
        curValues = curValues(~isnan(curValues));

        nValues(i) = numel(curValues);

        if nValues(i) >= 1
            meanValue(i) = mean(curValues);
        end

        if nValues(i) >= 2
            stdValue(i) = std(curValues, 0);
            varValue(i) = var(curValues, 0);
            semValue(i) = stdValue(i) / sqrt(nValues(i));
        else
            % A division with one value has a mean but no estimable variance.
            stdValue(i) = NaN;
            varValue(i) = NaN;
            semValue(i) = NaN;
        end
    end

    % ------------------------------------------------------------
    % beta0 and alpha_i.
    % ------------------------------------------------------------
    refIdx = divisions == referenceDivision;

    beta0 = meanValue(refIdx);

    alphaEstimate = meanValue - beta0;
    percentBias = 100 * alphaEstimate ./ beta0;

    % ------------------------------------------------------------
    % Optional Welch/Satterthwaite uncertainty for alpha_i.
    %
    % This requires variance estimates in both groups.
    % If referenceDivision = 36, n = 1, so pAlpha will usually be NaN.
    % That is expected and statistically honest.
    % ------------------------------------------------------------
    refVar = varValue(refIdx);
    refN = nValues(refIdx);

    seAlpha = nan(nDivisions, 1);
    dfAlpha = nan(nDivisions, 1);
    tAlpha = nan(nDivisions, 1);
    pAlpha = nan(nDivisions, 1);

    for i = 1:nDivisions

        if divisions(i) == referenceDivision
            seAlpha(i) = 0;
            continue;
        end

        curVar = varValue(i);
        curN = nValues(i);

        if curN >= 2 && refN >= 2 && ~isnan(curVar) && ~isnan(refVar)

            a = curVar / curN;
            b = refVar / refN;

            seAlpha(i) = sqrt(a + b);

            dfAlpha(i) = (a + b)^2 / ...
                (a^2 / (curN - 1) + b^2 / (refN - 1));

            tAlpha(i) = alphaEstimate(i) / seAlpha(i);
            pAlpha(i) = 2 * tcdf(-abs(tAlpha(i)), dfAlpha(i));
        end
    end

    % ------------------------------------------------------------
    % Output table.
    % ------------------------------------------------------------
    biasTable = table( ...
        divisions, ...
        nValues, ...
        meanValue, ...
        stdValue, ...
        varValue, ...
        semValue, ...
        repmat(beta0, nDivisions, 1), ...
        alphaEstimate, ...
        percentBias, ...
        seAlpha, ...
        dfAlpha, ...
        tAlpha, ...
        pAlpha, ...
        'VariableNames', { ...
        'division', ...
        'nValues', ...
        'meanValue', ...
        'stdValue', ...
        'varValue', ...
        'semValue', ...
        'beta0ReferenceMean', ...
        'alphaEstimate', ...
        'percentBias', ...
        'seAlpha', ...
        'dfAlpha', ...
        'tAlpha', ...
        'pAlpha'});

    result = struct;
    result.metabolite = metabName;
    result.modelText = "X_ij = beta0 + alpha_i + epsilon_ij, epsilon_ij ~ N(0, sigma_i^2)";
    result.referenceDivision = referenceDivision;
    result.beta0 = beta0;
    result.biasTable = biasTable;
    result.data = Tm;
end

function PlotDivisionHeteroscedasticBias(result)
% PlotDivisionHeteroscedasticBias
%
% Plots alpha_i as percent bias relative to the reference division.

    T = result.biasTable;

    x = 1:height(T);
    xLabels = string(T.division);

    figure;
    hold on;

    yline(0, '--', ...
        'LineWidth', 1.2, ...
        'DisplayName', 'Reference / no bias');

    plot(x, T.percentBias, '-o', ...
        'LineWidth', 2, ...
        'MarkerSize', 7, ...
        'DisplayName', 'Estimated fixed-effect bias');

    xticks(x);
    xticklabels(xLabels);

    xlabel('Scans per group / division size');
    ylabel(sprintf('Percent bias relative to Division_%d', ...
        result.referenceDivision));

    title(sprintf('%s division fixed-effect bias', result.metabolite), ...
        'Interpreter', 'none');

    legend('Location', 'best', 'Interpreter', 'none');
    grid on;
    box on;

    hold off;
end


function coordDir = ResolvePatientCoordDir(baseCoordDir, patientID, usePatientSubfolders)
% ResolvePatientCoordDir
%
% If usePatientSubfolders is true, returns:
%   fullfile(baseCoordDir, safePatientID)
%
% But if baseCoordDir already ends in the patient ID, it is left unchanged.
% This prevents accidentally reading from LCMFit/P01/P01.

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

    runOptions.inputMode = "singleFile";
    runOptions.dataDir = "";
    runOptions.twixFile = "";
    runOptions.filePattern = "*.dat";
    runOptions.recursive = false;
    runOptions.mustContain = strings(0, 1);
    runOptions.mustNotContain = strings(0, 1);
    runOptions.patientIDMode = "sequential";
    runOptions.singlePatientID = "P01";
    runOptions.coordDir = "";
    runOptions.loadMode = "allSubfolders";
    runOptions.selectedPatientID = "";
    runOptions.selectedSubfolder = "";

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

    if HasFieldOrProperty(cfg, "load")
        cfgLoad = GetFieldOrProperty(cfg, "load");

        if HasFieldOrProperty(cfgLoad, "mode")
            runOptions.loadMode = string(GetFieldOrProperty(cfgLoad, "mode"));
        end

        if HasFieldOrProperty(cfgLoad, "coordDir")
            runOptions.coordDir = string(GetFieldOrProperty(cfgLoad, "coordDir"));
        end

        if HasFieldOrProperty(cfgLoad, "selectedPatientID")
            runOptions.selectedPatientID = string(GetFieldOrProperty(cfgLoad, "selectedPatientID"));
        end

        if HasFieldOrProperty(cfgLoad, "selectedSubfolder")
            runOptions.selectedSubfolder = string(GetFieldOrProperty(cfgLoad, "selectedSubfolder"));
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
    
    p.addParameter('loadMode', runOptions.loadMode, @(x) ischar(x) || isstring(x));
    p.addParameter('selectedPatientID', runOptions.selectedPatientID, @(x) ischar(x) || isstring(x));
    p.addParameter('selectedSubfolder', runOptions.selectedSubfolder, @(x) ischar(x) || isstring(x));

    p.addParameter('usePatientSubfolders', runOptions.usePatientSubfolders, @(x) islogical(x) && isscalar(x));    p.addParameter('addPatientPrefixToFilenames', runOptions.addPatientPrefixToFilenames, @(x) islogical(x) && isscalar(x));
    p.addParameter('deleteOldLCModelFiles', runOptions.deleteOldLCModelFiles, @(x) islogical(x) && isscalar(x));
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

    runOptions.loadMode = string(p.Results.loadMode);
    runOptions.selectedPatientID = string(p.Results.selectedPatientID);
    runOptions.selectedSubfolder = string(p.Results.selectedSubfolder);

    runOptions.usePatientSubfolders = logical(p.Results.usePatientSubfolders);
    runOptions.addPatientPrefixToFilenames = logical(p.Results.addPatientPrefixToFilenames);
    runOptions.deleteOldLCModelFiles = logical(p.Results.deleteOldLCModelFiles);
end


function stats = InitializeCoordCacheStats(cacheOptions)

    stats = struct();
    stats.enabled = cacheOptions.enabled;
    stats.hits = 0;
    stats.missing = 0;
    stats.invalidated = 0;
    stats.corrupt = 0;
    stats.forced = 0;
    stats.disabled = 0;
    stats.writeFailures = 0;
    stats.validationSeconds = 0;
    stats.loadSeconds = 0;
    stats.parseSeconds = 0;
    stats.writeSeconds = 0;
    stats.coordStageSeconds = 0;
    stats.parserCallCount = 0;
end


function stats = UpdateCoordCacheStats(stats, event)

    switch event.category
        case "hit"
            stats.hits = stats.hits + 1;
        case "missing"
            stats.missing = stats.missing + 1;
        case "invalidated"
            stats.invalidated = stats.invalidated + 1;
        case "corrupt"
            stats.corrupt = stats.corrupt + 1;
        case "forced"
            stats.forced = stats.forced + 1;
        case "disabled"
            stats.disabled = stats.disabled + 1;
    end

    if event.didAttemptWrite && ~event.writeSucceeded
        stats.writeFailures = stats.writeFailures + 1;
    end

    stats.validationSeconds = stats.validationSeconds + event.validationSeconds;
    stats.loadSeconds = stats.loadSeconds + event.loadSeconds;
    stats.parseSeconds = stats.parseSeconds + event.parseSeconds;
    stats.writeSeconds = stats.writeSeconds + event.writeSeconds;
    stats.coordStageSeconds = stats.coordStageSeconds + event.coordStageSeconds;
    stats.parserCallCount = stats.parserCallCount + event.parserCallCount;
end


function PrintCoordCacheSummary(stats, totalSeconds)

    if stats.enabled
        fprintf(['\nCoord cache: %d hits, %d misses, %d invalidated, ', ...
            '%d corrupt, %d forced refresh, %d write failures.\n'], ...
            stats.hits, ...
            stats.missing, ...
            stats.invalidated, ...
            stats.corrupt, ...
            stats.forced, ...
            stats.writeFailures);
    else
        fprintf('\nCoord cache: disabled; %d patient(s) parsed normally.\n', stats.disabled);
    end

    downstreamSeconds = max(totalSeconds - stats.coordStageSeconds, 0);
    fprintf(['Coord timing (seconds): validation %.3f, cache load %.3f, ', ...
        'parse %.3f, cache write %.3f, downstream %.3f, total %.3f.\n'], ...
        stats.validationSeconds, ...
        stats.loadSeconds, ...
        stats.parseSeconds, ...
        stats.writeSeconds, ...
        downstreamSeconds, ...
        totalSeconds);
    fprintf('VDIIO.ReadLCMCoord calls: %d.\n', stats.parserCallCount);
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

function tablesByMetric = BuildTablesByMetric(coordTable, metabList, metrics)
% BuildTablesByMetric
%
% Creates:
%   tablesByMetric.sig
%   tablesByMetric.ratioCr
%   tablesByMetric.CRLB

    metrics = string(metrics(:));
    tablesByMetric = struct;

    for i = 1:numel(metrics)

        metricName = metrics(i);
        metricField = matlab.lang.makeValidName(metricName);

        tablesByMetric.(metricField) = coordToMetabs( ...
            coordTable, ...
            metabList, ...
            metricName);
    end
end

function [allResults, skippedTable] = FitAllMetabolitesSinglePatient(modelTable, metabList, varargin)
% FitAllMetabolitesSinglePatient
%
% Runs the single-patient heteroscedastic fixed-effect model for every
% metabolite. Metabolites with no valid reference value are skipped.

    p = inputParser;
    p.addParameter('referenceDivision', 36, @(x) isnumeric(x) && isscalar(x));
    parse(p, varargin{:});

    referenceDivision = p.Results.referenceDivision;
    metabList = string(metabList(:));

    allResults = struct;

    skippedMetabolite = strings(0, 1);
    skippedReason = strings(0, 1);

    for i = 1:numel(metabList)

        metabName = metabList(i);
        fieldName = matlab.lang.makeValidName(metabName);

        try
            allResults.(fieldName) = FitDivisionHeteroscedasticFixedEffect( ...
                modelTable, ...
                metabName, ...
                'referenceDivision', referenceDivision);

        catch ME
            skippedMetabolite(end+1, 1) = metabName; %#ok<AGROW>
            skippedReason(end+1, 1) = string(ME.message); %#ok<AGROW>
        end
    end

    skippedTable = table(skippedMetabolite, skippedReason, ...
        'VariableNames', {'metabolite', 'reason'});
end

function metricTables = GetMetricTablesFromPatientResults(patientResults, metricName)
% GetMetricTablesFromPatientResults
%
% Returns a cell array with one tablesByDivision struct per patient for the
% requested metric.

    metricName = string(metricName);
    metricField = matlab.lang.makeValidName(metricName);

    nPatients = numel(patientResults);
    metricTables = cell(nPatients, 1);

    for pIdx = 1:nPatients

        if ~isfield(patientResults(pIdx).tablesByMetric, metricField)
            error('Patient %d is missing metric field "%s".', pIdx, metricField);
        end

        metricTables{pIdx} = patientResults(pIdx).tablesByMetric.(metricField);
    end
end

function [patientBiasTable, skippedTable] = ComputeMultiPatientDivisionBias(patientTablesByDivision, metabList, varargin)
% ComputeMultiPatientDivisionBias
%
% Computes division bias inside each patient:
%
%   alpha_p,i = mean(patient p, Division_i) - mean(patient p, reference division)
%
% Output has one row per:
%   patient x metabolite x division

    p = inputParser;
    p.addParameter('patientIDs', strings(0, 1), @(x) isstring(x) || iscellstr(x) || ischar(x));
    p.addParameter('referenceDivision', 36, @(x) isnumeric(x) && isscalar(x));
    p.addParameter('ignoreZeros', true, @(x) islogical(x) && isscalar(x));
    p.addParameter('metricName', "sig", @(x) ischar(x) || isstring(x));
    parse(p, varargin{:});

    patientIDs = string(p.Results.patientIDs);
    referenceDivision = p.Results.referenceDivision;
    ignoreZeros = p.Results.ignoreZeros;
    metricName = string(p.Results.metricName);

    metabList = string(metabList(:));

    if ~iscell(patientTablesByDivision)
        error('patientTablesByDivision must be a cell array, one entry per patient.');
    end

    nPatients = numel(patientTablesByDivision);

    if isempty(patientIDs)
        patientIDs = "P" + compose("%02d", (1:nPatients)');
    end

    patientIDs = patientIDs(:);

    if numel(patientIDs) ~= nPatients
        error('Number of patientIDs must match number of patient tables.');
    end

    patientIndexCol = zeros(0, 1);
    patientIDCol = strings(0, 1);
    metricCol = strings(0, 1);
    metaboliteCol = strings(0, 1);
    divisionCol = zeros(0, 1);
    nValuesCol = zeros(0, 1);
    meanValueCol = zeros(0, 1);
    stdValueCol = zeros(0, 1);
    varValueCol = zeros(0, 1);
    semValueCol = zeros(0, 1);
    beta0Col = zeros(0, 1);
    alphaCol = zeros(0, 1);
    percentCol = zeros(0, 1);

    skippedPatientIndex = zeros(0, 1);
    skippedPatientID = strings(0, 1);
    skippedMetric = strings(0, 1);
    skippedMetabolite = strings(0, 1);
    skippedReason = strings(0, 1);

    for pIdx = 1:nPatients

        curPatientTables = patientTablesByDivision{pIdx};

        modelTable = BuildConcentrationModelTable( ...
            curPatientTables, ...
            metabList, ...
            'ignoreZeros', ignoreZeros);

        for m = 1:numel(metabList)

            metabName = metabList(m);

            try
                curResult = FitDivisionHeteroscedasticFixedEffect( ...
                    modelTable, ...
                    metabName, ...
                    'referenceDivision', referenceDivision);

                T = curResult.biasTable;
                nRows = height(T);

                patientIndexCol = [patientIndexCol; repmat(pIdx, nRows, 1)]; %#ok<AGROW>
                patientIDCol = [patientIDCol; repmat(patientIDs(pIdx), nRows, 1)]; %#ok<AGROW>
                metricCol = [metricCol; repmat(metricName, nRows, 1)]; %#ok<AGROW>
                metaboliteCol = [metaboliteCol; repmat(metabName, nRows, 1)]; %#ok<AGROW>
                divisionCol = [divisionCol; T.division]; %#ok<AGROW>
                nValuesCol = [nValuesCol; T.nValues]; %#ok<AGROW>
                meanValueCol = [meanValueCol; T.meanValue]; %#ok<AGROW>
                stdValueCol = [stdValueCol; T.stdValue]; %#ok<AGROW>
                varValueCol = [varValueCol; T.varValue]; %#ok<AGROW>
                semValueCol = [semValueCol; T.semValue]; %#ok<AGROW>
                beta0Col = [beta0Col; T.beta0ReferenceMean]; %#ok<AGROW>
                alphaCol = [alphaCol; T.alphaEstimate]; %#ok<AGROW>
                percentCol = [percentCol; T.percentBias]; %#ok<AGROW>

            catch ME
                skippedPatientIndex(end+1, 1) = pIdx; %#ok<AGROW>
                skippedPatientID(end+1, 1) = patientIDs(pIdx); %#ok<AGROW>
                skippedMetric(end+1, 1) = metricName; %#ok<AGROW>
                skippedMetabolite(end+1, 1) = metabName; %#ok<AGROW>
                skippedReason(end+1, 1) = string(ME.message); %#ok<AGROW>
            end
        end
    end

    patientBiasTable = table( ...
        patientIndexCol, ...
        patientIDCol, ...
        metricCol, ...
        metaboliteCol, ...
        divisionCol, ...
        nValuesCol, ...
        meanValueCol, ...
        stdValueCol, ...
        varValueCol, ...
        semValueCol, ...
        beta0Col, ...
        alphaCol, ...
        percentCol, ...
        'VariableNames', { ...
        'patientIndex', ...
        'patientID', ...
        'metric', ...
        'metabolite', ...
        'division', ...
        'nValues', ...
        'meanValue', ...
        'stdValue', ...
        'varValue', ...
        'semValue', ...
        'beta0ReferenceMean', ...
        'alphaEstimate', ...
        'percentBias'});

    skippedTable = table( ...
        skippedPatientIndex, ...
        skippedPatientID, ...
        skippedMetric, ...
        skippedMetabolite, ...
        skippedReason, ...
        'VariableNames', { ...
        'patientIndex', ...
        'patientID', ...
        'metric', ...
        'metabolite', ...
        'reason'});
end

function summaryTable = SummarizeMultiPatientDivisionBias(patientBiasTable, varargin)
% SummarizeMultiPatientDivisionBias
%
% Summarizes alpha_i and percent bias across patients.
%
% For each metabolite and division, this tests:
%   H0: mean patient-level bias = 0

    p = inputParser;
    p.addParameter('referenceDivision', 36, @(x) isnumeric(x) && isscalar(x));
    p.addParameter('alphaLevel', 0.05, @(x) isnumeric(x) && isscalar(x));
    parse(p, varargin{:});

    referenceDivision = p.Results.referenceDivision;
    alphaLevel = p.Results.alphaLevel;

    if isempty(patientBiasTable)
        summaryTable = MakeEmptySummaryBiasTable();
        return;
    end

    metrics = unique(patientBiasTable.metric, 'stable');
    metabolites = unique(patientBiasTable.metabolite, 'stable');
    divisions = unique(patientBiasTable.division);
    divisions = sort(divisions(:));

    metricCol = strings(0, 1);
    metaboliteCol = strings(0, 1);
    divisionCol = zeros(0, 1);

    nPatientsAlphaCol = zeros(0, 1);
    meanAlphaCol = zeros(0, 1);
    stdAlphaCol = zeros(0, 1);
    semAlphaCol = zeros(0, 1);
    tAlphaCol = zeros(0, 1);
    dfAlphaCol = zeros(0, 1);
    pAlphaCol = zeros(0, 1);
    alphaCILowCol = zeros(0, 1);
    alphaCIHighCol = zeros(0, 1);

    nPatientsPercentCol = zeros(0, 1);
    meanPercentCol = zeros(0, 1);
    stdPercentCol = zeros(0, 1);
    semPercentCol = zeros(0, 1);
    tPercentCol = zeros(0, 1);
    dfPercentCol = zeros(0, 1);
    pPercentCol = zeros(0, 1);
    percentCILowCol = zeros(0, 1);
    percentCIHighCol = zeros(0, 1);

    isReferenceCol = false(0, 1);

    for metricIdx = 1:numel(metrics)

        metricName = metrics(metricIdx);

        for m = 1:numel(metabolites)

            metabName = metabolites(m);

            for d = 1:numel(divisions)

                curDivision = divisions(d);

                idx = patientBiasTable.metric == metricName & ...
                      patientBiasTable.metabolite == metabName & ...
                      patientBiasTable.division == curDivision;

                curRows = patientBiasTable(idx, :);

                alphaVals = curRows.alphaEstimate;
                percentVals = curRows.percentBias;

                [meanAlpha, stdAlpha, semAlpha, tAlpha, dfAlpha, pAlpha, ciAlphaLow, ciAlphaHigh, nAlpha] = ...
                    OneSampleSummary(alphaVals, alphaLevel);

                [meanPercent, stdPercent, semPercent, tPercent, dfPercent, pPercent, ciPercentLow, ciPercentHigh, nPercent] = ...
                    OneSampleSummary(percentVals, alphaLevel);

                metricCol(end+1, 1) = metricName; %#ok<AGROW>
                metaboliteCol(end+1, 1) = metabName; %#ok<AGROW>
                divisionCol(end+1, 1) = curDivision; %#ok<AGROW>

                nPatientsAlphaCol(end+1, 1) = nAlpha; %#ok<AGROW>
                meanAlphaCol(end+1, 1) = meanAlpha; %#ok<AGROW>
                stdAlphaCol(end+1, 1) = stdAlpha; %#ok<AGROW>
                semAlphaCol(end+1, 1) = semAlpha; %#ok<AGROW>
                tAlphaCol(end+1, 1) = tAlpha; %#ok<AGROW>
                dfAlphaCol(end+1, 1) = dfAlpha; %#ok<AGROW>
                pAlphaCol(end+1, 1) = pAlpha; %#ok<AGROW>
                alphaCILowCol(end+1, 1) = ciAlphaLow; %#ok<AGROW>
                alphaCIHighCol(end+1, 1) = ciAlphaHigh; %#ok<AGROW>

                nPatientsPercentCol(end+1, 1) = nPercent; %#ok<AGROW>
                meanPercentCol(end+1, 1) = meanPercent; %#ok<AGROW>
                stdPercentCol(end+1, 1) = stdPercent; %#ok<AGROW>
                semPercentCol(end+1, 1) = semPercent; %#ok<AGROW>
                tPercentCol(end+1, 1) = tPercent; %#ok<AGROW>
                dfPercentCol(end+1, 1) = dfPercent; %#ok<AGROW>
                pPercentCol(end+1, 1) = pPercent; %#ok<AGROW>
                percentCILowCol(end+1, 1) = ciPercentLow; %#ok<AGROW>
                percentCIHighCol(end+1, 1) = ciPercentHigh; %#ok<AGROW>

                isReferenceCol(end+1, 1) = curDivision == referenceDivision; %#ok<AGROW>
            end
        end
    end

    summaryTable = table( ...
        metricCol, ...
        metaboliteCol, ...
        divisionCol, ...
        nPatientsAlphaCol, ...
        meanAlphaCol, ...
        stdAlphaCol, ...
        semAlphaCol, ...
        tAlphaCol, ...
        dfAlphaCol, ...
        pAlphaCol, ...
        alphaCILowCol, ...
        alphaCIHighCol, ...
        nPatientsPercentCol, ...
        meanPercentCol, ...
        stdPercentCol, ...
        semPercentCol, ...
        tPercentCol, ...
        dfPercentCol, ...
        pPercentCol, ...
        percentCILowCol, ...
        percentCIHighCol, ...
        isReferenceCol, ...
        'VariableNames', { ...
        'metric', ...
        'metabolite', ...
        'division', ...
        'nPatientsAlpha', ...
        'meanAlpha', ...
        'stdAlpha', ...
        'semAlpha', ...
        'tAlpha', ...
        'dfAlpha', ...
        'pAlpha', ...
        'alphaCILow', ...
        'alphaCIHigh', ...
        'nPatientsPercent', ...
        'meanPercentBias', ...
        'stdPercentBias', ...
        'semPercentBias', ...
        'tPercentBias', ...
        'dfPercentBias', ...
        'pPercentBias', ...
        'percentCILow', ...
        'percentCIHigh', ...
        'isReferenceDivision'});
end

function summaryTable = MakeEmptySummaryBiasTable()

    summaryTable = table( ...
        strings(0, 1), strings(0, 1), zeros(0, 1), ...
        zeros(0, 1), zeros(0, 1), zeros(0, 1), zeros(0, 1), ...
        zeros(0, 1), zeros(0, 1), zeros(0, 1), zeros(0, 1), zeros(0, 1), ...
        zeros(0, 1), zeros(0, 1), zeros(0, 1), zeros(0, 1), ...
        zeros(0, 1), zeros(0, 1), zeros(0, 1), zeros(0, 1), zeros(0, 1), ...
        false(0, 1), ...
        'VariableNames', { ...
        'metric', 'metabolite', 'division', ...
        'nPatientsAlpha', 'meanAlpha', 'stdAlpha', 'semAlpha', ...
        'tAlpha', 'dfAlpha', 'pAlpha', 'alphaCILow', 'alphaCIHigh', ...
        'nPatientsPercent', 'meanPercentBias', 'stdPercentBias', 'semPercentBias', ...
        'tPercentBias', 'dfPercentBias', 'pPercentBias', 'percentCILow', 'percentCIHigh', ...
        'isReferenceDivision'});
end

function [meanVal, stdVal, semVal, tVal, dfVal, pVal, ciLow, ciHigh, n] = OneSampleSummary(values, alphaLevel)
% OneSampleSummary
%
% One-sample t summary for testing mean(values) = 0.
% The independent unit here should be patient-level bias estimates.

    values = values(:);
    values = values(~isnan(values));

    n = numel(values);

    meanVal = NaN;
    stdVal = NaN;
    semVal = NaN;
    tVal = NaN;
    dfVal = NaN;
    pVal = NaN;
    ciLow = NaN;
    ciHigh = NaN;

    if n == 0
        return;
    end

    meanVal = mean(values);

    if n == 1
        return;
    end

    stdVal = std(values, 0);
    semVal = stdVal / sqrt(n);
    dfVal = n - 1;

    if semVal == 0
        if meanVal == 0
            tVal = 0;
            pVal = 1;
        else
            tVal = sign(meanVal) * Inf;
            pVal = 0;
        end

        ciLow = meanVal;
        ciHigh = meanVal;
        return;
    end

    tVal = meanVal / semVal;
    pVal = 2 * tcdf(-abs(tVal), dfVal);

    tCrit = tinv(1 - alphaLevel/2, dfVal);

    ciLow = meanVal - tCrit * semVal;
    ciHigh = meanVal + tCrit * semVal;
end

function PlotMultiPatientBiasSummary(summaryTable, metabName, varargin)
% PlotMultiPatientBiasSummary
%
% Plots mean percent bias across patients.

    p = inputParser;
    p.addParameter('metricName', "sig", @(x) ischar(x) || isstring(x));
    p.addParameter('excludeReference', false, @(x) islogical(x) && isscalar(x));
    parse(p, varargin{:});

    metricName = string(p.Results.metricName);
    excludeReference = p.Results.excludeReference;

    metabName = string(metabName);

    T = summaryTable(strcmpi(summaryTable.metabolite, metabName) & ...
                     summaryTable.metric == metricName, :);

    if excludeReference
        T = T(~T.isReferenceDivision, :);
    end

    if isempty(T)
        error('No rows found for metabolite "%s" and metric "%s".', metabName, metricName);
    end

    T = sortrows(T, 'division');

    x = 1:height(T);

    figure;
    hold on;

    yline(0, '--', ...
        'LineWidth', 1.2, ...
        'DisplayName', 'No bias');

    errorbar(x, T.meanPercentBias, T.semPercentBias, '-o', ...
        'LineWidth', 2, ...
        'MarkerSize', 7, ...
        'DisplayName', 'Mean percent bias +/- SEM');

    xticks(x);
    xticklabels(string(T.division));

    xlabel('Scans per group / division size');
    ylabel('Mean percent bias across patients');

    title(sprintf('%s multi-patient division bias (%s)', metabName, metricName), ...
        'Interpreter', 'none');

    legend('Location', 'best', 'Interpreter', 'none');
    grid on;
    box on;

    hold off;
end

function PlotMultiPatientBiasSpaghetti(patientBiasTable, metabName, varargin)
% PlotMultiPatientBiasSpaghetti
%
% One line per patient, showing patient-level percent bias across divisions.

    p = inputParser;
    p.addParameter('metricName', "sig", @(x) ischar(x) || isstring(x));
    p.addParameter('excludeReference', false, @(x) islogical(x) && isscalar(x));
    parse(p, varargin{:});

    metricName = string(p.Results.metricName);
    excludeReference = p.Results.excludeReference;
    metabName = string(metabName);

    T = patientBiasTable(strcmpi(patientBiasTable.metabolite, metabName) & ...
                         patientBiasTable.metric == metricName, :);

    if excludeReference
        T = T(T.percentBias ~= 0 | T.division ~= max(T.division), :);
    end

    if isempty(T)
        error('No rows found for metabolite "%s" and metric "%s".', metabName, metricName);
    end

    divisions = sort(unique(T.division));
    patientIDs = unique(T.patientID, 'stable');

    figure;
    hold on;

    for pIdx = 1:numel(patientIDs)

        curPatient = patientIDs(pIdx);
        y = nan(numel(divisions), 1);

        for d = 1:numel(divisions)
            idx = T.patientID == curPatient & T.division == divisions(d);
            if any(idx)
                y(d) = T.percentBias(find(idx, 1, 'first'));
            end
        end

        plot(1:numel(divisions), y, '-o', ...
            'LineWidth', 1.2, ...
            'DisplayName', curPatient);
    end

    yline(0, '--', ...
        'LineWidth', 1.2, ...
        'DisplayName', 'No bias');

    xticks(1:numel(divisions));
    xticklabels(string(divisions));

    xlabel('Scans per group / division size');
    ylabel('Percent bias relative to patient reference');

    title(sprintf('%s patient-level division bias (%s)', metabName, metricName), ...
        'Interpreter', 'none');

    legend('Location', 'best', 'Interpreter', 'none');
    grid on;
    box on;

    hold off;
end

function [covTable, corrTable, biasMatrixTable] = ComputeBiasCovarianceForMetabolite(patientBiasTable, metabName, varargin)
% ComputeBiasCovarianceForMetabolite
%
% Builds a patient x division matrix and computes pairwise covariance and
% correlation across divisions.
%
% Rows of the bias matrix are patients.
% Columns are divisions.

    p = inputParser;
    p.addParameter('metricName', "sig", @(x) ischar(x) || isstring(x));
    p.addParameter('valueColumn', "percentBias", @(x) ischar(x) || isstring(x));
    p.addParameter('excludeReference', false, @(x) islogical(x) && isscalar(x));
    parse(p, varargin{:});

    metricName = string(p.Results.metricName);
    valueColumn = string(p.Results.valueColumn);
    excludeReference = p.Results.excludeReference;

    if ~ismember(valueColumn, string(patientBiasTable.Properties.VariableNames))
        error('patientBiasTable does not contain value column "%s".', valueColumn);
    end

    metabName = string(metabName);

    T = patientBiasTable(strcmpi(patientBiasTable.metabolite, metabName) & ...
                         patientBiasTable.metric == metricName, :);

    if excludeReference
        T = T(T.(valueColumn) ~= 0 | T.division ~= max(T.division), :);
    end

    if isempty(T)
        error('No rows found for metabolite "%s" and metric "%s".', metabName, metricName);
    end

    patientIDs = unique(T.patientID, 'stable');
    divisions = sort(unique(T.division));

    A = nan(numel(patientIDs), numel(divisions));

    for pIdx = 1:numel(patientIDs)
        for d = 1:numel(divisions)
            idx = T.patientID == patientIDs(pIdx) & T.division == divisions(d);
            if any(idx)
                vals = T.(valueColumn)(idx);
                A(pIdx, d) = vals(1);
            end
        end
    end

    C = PairwiseCovariance(A);
    R = PairwiseCorrelation(A);

    varNames = matlab.lang.makeValidName("Division_" + string(divisions));

    covTable = array2table(C, ...
        'VariableNames', cellstr(varNames), ...
        'RowNames', cellstr(varNames));

    corrTable = array2table(R, ...
        'VariableNames', cellstr(varNames), ...
        'RowNames', cellstr(varNames));

    biasMatrixTable = array2table(A, ...
        'VariableNames', cellstr(varNames), ...
        'RowNames', cellstr(patientIDs));
end

function C = PairwiseCovariance(A)

    nCols = size(A, 2);
    C = nan(nCols, nCols);

    for i = 1:nCols
        for j = 1:nCols
            keep = ~isnan(A(:, i)) & ~isnan(A(:, j));
            if sum(keep) >= 2
                tmp = cov(A(keep, i), A(keep, j));
                C(i, j) = tmp(1, 2);
            end
        end
    end
end

function R = PairwiseCorrelation(A)

    nCols = size(A, 2);
    R = nan(nCols, nCols);

    for i = 1:nCols
        for j = 1:nCols
            keep = ~isnan(A(:, i)) & ~isnan(A(:, j));
            if sum(keep) >= 2
                tmp = corrcoef(A(keep, i), A(keep, j));
                R(i, j) = tmp(1, 2);
            end
        end
    end
end

function patients = BuildPatientTableForLoading(cfg, runOptions)

loadMode = lower(string(runOptions.loadMode));
baseCoordDir = string(runOptions.coordDir);

if strlength(baseCoordDir) == 0
    error('coordDir is empty.');
end

if ~isfolder(baseCoordDir)
    error('coordDir does not exist: %s', baseCoordDir);
end

switch loadMode

    case "allsubfolders"

        info = dir(baseCoordDir);
        info = info([info.isdir]);

        folderNames = string({info.name});
        keep = folderNames ~= "." & folderNames ~= "..";
        info = info(keep);

        patientID = strings(0, 1);
        twixFile = strings(0, 1);
        coordDir = strings(0, 1);

        for k = 1:numel(info)

            curFolder = string(fullfile(info(k).folder, info(k).name));
            coordInfo = dir(fullfile(curFolder, "*.coord"));

            if isempty(coordInfo)
                continue;
            end

            curPatientID = string(matlab.lang.makeValidName(info(k).name));

            patientID(end+1, 1) = curPatientID; %#ok<AGROW>
            twixFile(end+1, 1) = ""; %#ok<AGROW>
            coordDir(end+1, 1) = curFolder; %#ok<AGROW>
        end

        if isempty(patientID)
            error('No patient subfolders with .coord files were found in: %s', baseCoordDir);
        end

        patients = table(patientID, twixFile, coordDir, ...
            'VariableNames', {'patientID', 'twixFile', 'coordDir'});

    case "singlesubfolder"

        selectedSubfolder = string(runOptions.selectedSubfolder);
        selectedPatientID = string(runOptions.selectedPatientID);

        if strlength(selectedSubfolder) > 0

            coordDir = selectedSubfolder;

            if ~isfolder(coordDir)
                error('selectedSubfolder does not exist: %s', coordDir);
            end

            if strlength(selectedPatientID) == 0
                [~, folderName] = fileparts(char(coordDir));
                selectedPatientID = string(folderName);
            end

        else

            if strlength(selectedPatientID) == 0
                error('singleSubfolder mode requires selectedPatientID or selectedSubfolder.');
            end

            safePatientID = string(matlab.lang.makeValidName(char(selectedPatientID)));
            coordDir = string(fullfile(baseCoordDir, safePatientID));

            if ~isfolder(coordDir)
                error('Patient folder does not exist: %s', coordDir);
            end
        end

        coordInfo = dir(fullfile(coordDir, "*.coord"));

        if isempty(coordInfo)
            error('No .coord files found in selected folder: %s', coordDir);
        end

        patientID = string(matlab.lang.makeValidName(char(selectedPatientID)));
        twixFile = "";
        patients = table(patientID, twixFile, coordDir, ...
            'VariableNames', {'patientID', 'twixFile', 'coordDir'});

    otherwise
        error('Unknown loadMode "%s". Use "allSubfolders" or "singleSubfolder".', runOptions.loadMode);
end
end
