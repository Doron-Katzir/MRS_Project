%% Configuration

clear;
img = VDIIO.LoadTwix(['C:\Users\doronkatzir1\Desktop\Thesis_Lab\Data\' ...
    'meas_MID00090_FID32072_eja_svs_slaser_TE_80_r0.dat'], ...
    'isICEChop', true);
setSizes = [1, 2, 4, 6, 9, 12, 18, 36];
basisMetabList = ["NAA", "NAAG", "Cr", "PCr", "GPC", "PCh", "Glu", "Gln", ...
        "GABA", "GSH", "Tau", "Asc", "Glc", "Ace", "mI", "sI", "Asp", ...
        "Lac"];
sumMetabList = ["GPC+PCh", "NAA+NAAG", "Cr+PCr", "Glu+Gln"];
analysisMetabList = [basisMetabList, sumMetabList];
coordDir = 'C:\Users\doronkatzir1\Desktop\Thesis_Lab\LCMFit';


%% Preprocessing

img.AddCoils;
img.ChopPts("numPts", 4096);
img.FT;
img.Phase; % This will automatically (hopefully!) phase your data
% img.Average("dim", "set");

%% Creating a simple basis set

% % We will simulate a .basis file with a set of basis functions (metabolites)
% % using the VDI routines with some explanations along the way.
% % Simplest case: "single molecule" at x=y=z=0, with an idealized simple
% % sequence (usually PRESS)
% % All of the simulations are carried out using the SpinsJ class.
% spins = SpinsJ('B0', img.B0);
% % Let's populate 18 different metabolite types, which will be stored in
% % the spins.metab property
% metabList = ["NAA", "NAAG", "Cr", "PCr", "GPC", "PCh", "Glu", "Gln", ...
%     "GABA", "GSH", "Tau", "Asc", "Glc", "Ace", "mI", "sI", "Asp", "Lac"];
% spins.AddMetab(metabList);
% % Now will create an idealized PRESS sequence.
% % Use the actual TE from the header data in your VDIImageND object,
% % which can be read using the GetTwixHeaderReport method.
% TE = VDIIO.GetTwixHeaderReport(img.metadata.hdr).TE;
% seq = Sequence.GetIdealSequence("PRESS", "TE", TE);
% [spinsOut, TT] = seq.Apply(spins, "isVerbose", true);
% % The output will be a 1x18 array of TransitionTable objects, each of
% % which will contain the "results" of the simulation for each of the 18
% % metabolites. For example, TT(4) is for PCr:
% % freqHz amp phaseDeg PPM3T PPM7T
% % _______ ____ __________ _______ ______
% %
% % -490.67 0.75 -0.011245 0.83862 3.0337
% % -222.9 0.5 -0.0051084 2.935 3.9321
% % Create a VDI basis-set object
% myBasis = VDIBasis(TT, 'B0', img.B0, ...
%     'numAcqPts', img.numSpecPts, ...
%     'dwellTime', img.dwellTime);
% myBasis.ExportBasisToLCModel("Doron.basis", 'TE', TE);

%% Subdividing a dataset

% % You will need to manually subdivide img.data and create individual
% % VDIImageND objects, which will then be fit (all using the same basis set).
% img1 = img.Copy;
% img1.data = img1.data(:,:,:,:,1:18);
% img1.FitLCModel("Doron.basis", 'isVerbose', true);
% % TODO: Look at time courses of metabolites at different "temporal resolutions"
% % 1x36
% % 2x18
% % 4x9
% % 6x6
% % 9x4
% % 18x2
% % 36x1
% % For each you can plot the same time course (over the same x-axis time range!)

%% Functions

function [myBasis, TE] = createBasis(inputImg ,metList)
    spins = SpinsJ('B0', inputImg.B0);
    spins.AddMetab(metList);
    TE = VDIIO.GetTwixHeaderReport(inputImg.metadata.hdr).TE;
    seq = Sequence.GetIdealSequence("PRESS", "TE", TE);
    [spinsOut, TT] = seq.Apply(spins, "isVerbose", true);
    myBasis = VDIBasis(TT, 'B0', inputImg.B0, ...
        'numAcqPts', inputImg.numSpecPts, ...
        'dwellTime', inputImg.dwellTime);
end 


function imgBlocks = FitSubset(inputImg, nAvg, basisName)
    % Get parameters.
    setString = 'set';
    dimTypes = inputImg.GetDimType;
    idxSetDim = find(strcmpi(dimTypes, setString), 1);
    totalSets = size(inputImg.data, idxSetDim);
    nBlocks = totalSets / nAvg;

    % Check if desired number of sets to average is an exact divisor of the
    % total number of sets.
    if mod(totalSets, nAvg) ~= 0
        error(['nAvg must evenly divide the total number of sets. ' ...
            'totalSets = %d, nAvg = %d.'], totalSets, nAvg);
    end
    
    % Create copies of the image according to desired number of sets.
    %   Preallocate cell array of sliced VDIImageND objects
    imgBlocks = cell(nBlocks, 1);
    setRanges = zeros(nBlocks, 2);

    for k = 1:nBlocks

        firstSet = (k - 1) * nAvg + 1;
        lastSet  = k * nAvg;
        setRanges(k, :) = [firstSet, lastSet];
        sliceArgs = repmat({':'}, 1, numel(dimTypes));
        sliceArgs{idxSetDim} = sprintf('%d:%d', firstSet, lastSet);

        % Slice original image into a copied subset
        imgBlock = inputImg.Slice(sliceArgs{:});
        
        % Average this subset over the set dimension FIRST
        imgBlock = imgBlock.Average("dim", setString);
        
        % Fit the averaged block
        fitFileName = "Division_" + nAvg + "_" + "part_" + k + ".basis";
        
        imgBlock.FitLCModel(basisName, ...
            'isVerbose', true, ...
            'outputFilename', fitFileName);
        
        imgBlocks{k} = imgBlock;
    end
end

function [quant, fitData] = readCoord(coordDir)
    coordInfo = dir(fullfile(coordDir, "*.coord"));
    
    if isempty(coordInfo)
        error('No .coord files found in: %s', coordDir);
    end
    
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
%% Run main

% Create and export basis function
[basisFunc, TE] = createBasis(img, basisMetabList);
basisName = "For_division_" + erase(img.name, " ") + ".basis";
basisFunc.ExportBasisToLCModel(basisName, 'TE', TE);

% Splice data to subsets
for size = setSizes
    splicedImg = FitSubset(img, size, basisName);
end

% Read coord files
[quant, fitData] = readCoord(coordDir);
coordTable = quant.metabTable;
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

% Plot one chosen metabolite
chosenMetab = "NAA+NAAG";

summaryChosen = PlotMetabAcrossDivisions( ...
    tablesByDivision, ...
    chosenMetab, ...
    'yLabel', "sig", ...
    'makeFigure', true);



