%% Configuration

cfg = ProjectConfig();

%% Main
% Read coord files
[quant, fitData, coordFiles] = readCoord(cfg.paths.coordDir);
coordTable = quant.metabTable;
analysisMetabList = cfg.metabolites.analysis;
tablesByDivision = coordToMetabs(coordTable, analysisMetabList, "sig");

% Plot
summaryByMetab = struct;

for i = 1:numel(analysisMetabList)

    metabName = analysisMetabList(i);
    fieldName = matlab.lang.makeValidName(metabName);

    summaryByMetab.(fieldName) = PlotMetabAcrossDivisions( ...
        tablesByDivision, ...
        metabName, ...
        'yLabel', "sig", ...
        'makeFigure', false);

end

% Plot for all chosen metabolites
tablesByMetric.sig = coordToMetabs(coordTable, analysisMetabList, "sig");
tablesByMetric.ratioCr = coordToMetabs(coordTable, analysisMetabList, "ratioCr");
tablesByMetric.CRLB = coordToMetabs(coordTable, analysisMetabList, "CRLB");
sumMetabList = cfg.metabolites.sum;

for metab = sumMetabList
    chosenMetab = metab;
    
    summaryChosen = PlotMetabAcrossDivisions( ...
        tablesByDivision, ...
        chosenMetab, ...
        'yLabel', "sig", ...
        'makeFigure', true);

    grandTables = PlotMetabGrandAverageBarsMultiMetric( ...
        tablesByMetric, ...
        metab, ...
        'barWidth', 0.35, ...
        'yPaddingFrac', 0.30);
end

% Plot spectra
PlotDivisionFittedSpectraStack(fitData, coordFiles, 6);

% ANOVA test
anovaSig = RunExploratoryDivisionAnova( ...
    tablesByDivision, ...
    analysisMetabList, ...
    'ignoreZeros', true, ...
    'showPlots', false);

anovaSig = sortrows(anovaSig, 'pValue');
disp(anovaSig)

%% Functions

function [quant, fitData, coordFiles] = readCoord(coordDir)

    coordInfo = dir(fullfile(coordDir, "*.coord"));
    
    if isempty(coordInfo)
        error('No .coord files found in: %s', coordDir);
    end

    % Sort filenames so part order is reproducible
    fileNames = string({coordInfo.name});
    [~, order] = sort(fileNames);
    coordInfo = coordInfo(order);
    
    coordFiles = strings(numel(coordInfo), 1);
    
    for k = 1:numel(coordInfo)
        coordFiles(k) = fullfile(coordInfo(k).folder, coordInfo(k).name);
    end

    [quant, fitData] = VDIIO.ReadLCMCoord(coordFiles);
end

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

    % Extract division size and part number from filename
    % Supports:
    %   Division_18_1.basis.coord
    %   Division_18_part_1.basis.coord
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

function anovaData = BuildDivisionLongTable(tablesByDivisionInput, metabList, varargin)
    % BuildDivisionLongTable
    %
    % Converts tablesByDivision into long-format data:
    %
    %   patientIndex | metabolite | division | part | value
    %
    % Supports:
    %   single patient:
    %       BuildDivisionLongTable(tablesByDivision, metabList)
    %
    %   multiple patients:
    %       BuildDivisionLongTable({tablesP1, tablesP2}, metabList)
    
    p = inputParser;
    p.addParameter('ignoreZeros', true, @(x) islogical(x) && isscalar(x));
    parse(p, varargin{:});
    
    ignoreZeros = p.Results.ignoreZeros;
    metabList = string(metabList(:));
    
    % Convert input to cell array of patients
    if iscell(tablesByDivisionInput)
        patientTables = tablesByDivisionInput(:);
    elseif isstruct(tablesByDivisionInput) && numel(tablesByDivisionInput) > 1
        patientTables = num2cell(tablesByDivisionInput);
    elseif isstruct(tablesByDivisionInput)
        patientTables = {tablesByDivisionInput};
    else
        error('Input must be a tablesByDivision struct, struct array, or cell array of structs.');
    end
    
    rows = {};
    
    for pIdx = 1:numel(patientTables)
    
        curStruct = patientTables{pIdx};
        divisionFields = string(fieldnames(curStruct));
    
        for d = 1:numel(divisionFields)
    
            curField = divisionFields(d);
            curDivision = ParseDivisionNumber(curField);
    
            if isnan(curDivision)
                continue;
            end
    
            curTable = curStruct.(curField);
    
            for m = 1:numel(metabList)
    
                metabName = metabList(m);
    
                y = ExtractMetabVectorFromTable(curTable, metabName);
    
                if isempty(y)
                    continue;
                end
    
                y = double(y(:));
    
                for partIdx = 1:numel(y)
    
                    curValue = y(partIdx);
    
                    if ignoreZeros && curValue == 0
                        curValue = NaN;
                    end
    
                    rows(end+1, :) = { ...
                        pIdx, ...
                        metabName, ...
                        curDivision, ...
                        partIdx, ...
                        curValue}; %#ok<AGROW>
                end
            end
        end
    end
    
    anovaData = cell2table(rows, ...
        'VariableNames', { ...
        'patientIndex', ...
        'metabolite', ...
        'division', ...
        'part', ...
        'value'});
    
    anovaData.metabolite = string(anovaData.metabolite);
end

function anovaResults = RunExploratoryDivisionAnova(tablesByDivision, metabList, varargin)
    % RunExploratoryDivisionAnova
    %
    % Exploratory one-way ANOVA:
    %   value ~ division
    %
    % This treats division parts as observations.
    % For one subject this is exploratory, not definitive inference.
    
    p = inputParser;
    p.addParameter('ignoreZeros', true, @(x) islogical(x) && isscalar(x));
    p.addParameter('showPlots', false, @(x) islogical(x) && isscalar(x));
    parse(p, varargin{:});
    
    ignoreZeros = p.Results.ignoreZeros;
    showPlots = p.Results.showPlots;
    
    metabList = string(metabList(:));
    
    anovaData = BuildDivisionLongTable(tablesByDivision, metabList, ...
        'ignoreZeros', ignoreZeros);
    
    rows = {};
    
    for m = 1:numel(metabList)
    
        metabName = metabList(m);
    
        Tm = anovaData(anovaData.metabolite == metabName, :);
        Tm = Tm(~isnan(Tm.value), :);
    
        if height(Tm) < 2 || numel(unique(Tm.division)) < 2
            rows(end+1, :) = {metabName, NaN, NaN, NaN, NaN}; %#ok<AGROW>
            continue;
        end
    
        if showPlots
            [pVal, tbl, stats] = anova1(Tm.value, categorical(Tm.division));
        else
            [pVal, tbl, stats] = anova1(Tm.value, categorical(Tm.division), 'off');
        end
    
        % Extract F statistic if available
        try
            F = tbl{2, 5};
            dfBetween = tbl{2, 3};
            dfWithin = tbl{3, 3};
        catch
            F = NaN;
            dfBetween = NaN;
            dfWithin = NaN;
        end
    
        rows(end+1, :) = { ...
            metabName, ...
            pVal, ...
            F, ...
            dfBetween, ...
            dfWithin}; %#ok<AGROW>
    end
    
    anovaResults = cell2table(rows, ...
        'VariableNames', { ...
        'metabolite', ...
        'pValue', ...
        'F', ...
        'dfBetween', ...
        'dfWithin'});
    
    anovaResults.metabolite = string(anovaResults.metabolite);
end