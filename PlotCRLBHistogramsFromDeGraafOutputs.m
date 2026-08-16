function crlbOutputs = PlotCRLBHistogramsFromDeGraafOutputs(deGraafOutputs, crlbCfg)
% PlotCRLBHistogramsFromDeGraafOutputs
%
% Collects raw CRLB values across patients and Division_1 parts, then plots
% one CRLB histogram per metabolite.
%
% Expected idea:
%   rows = patient/part/metabolite observations
%   columns = patientID, part, metabolite, CRLB

if nargin < 2 || isempty(crlbCfg)
    crlbCfg = struct();
end

opts = ParseCRLBHistogramOptions(crlbCfg);

if ~isfield(deGraafOutputs, "patientResultsByID")
    error("deGraafOutputs.patientResultsByID is missing.");
end

patientFields = string(fieldnames(deGraafOutputs.patientResultsByID));

crlbLongTable = table();

for pIdx = 1:numel(patientFields)

    patientField = patientFields(pIdx);
    patientResult = deGraafOutputs.patientResultsByID.(patientField);

    patientTable = ExtractPatientCRLBLongTable(patientResult, patientField);

    crlbLongTable = [crlbLongTable; patientTable]; %#ok<AGROW>
end

if isempty(crlbLongTable)
    error("No CRLB values were found. The DeGraaf output may not currently store raw part-level CRLB values.");
end

% Select metabolites.
if isscalar(string(opts.metabolites)) && strcmpi(string(opts.metabolites), "all")
    selectedMetabs = unique(string(crlbLongTable.metabolite), "stable");
else
    selectedMetabs = string(opts.metabolites(:));
end

% Apply sum-preferred filtering if requested.
if opts.useSumPreferredMetabolites
    selectedMetabs = ApplySumPreferredFiltering(selectedMetabs, opts.sumMetabolites);
end

% Keep only selected metabolites.
keepRows = ismember(string(crlbLongTable.metabolite), selectedMetabs);
crlbLongTable = crlbLongTable(keepRows, :);

if isempty(crlbLongTable)
    error("No CRLB values remained after metabolite filtering.");
end

% Create output directory.
if opts.saveFigures && ~exist(opts.outputDir, "dir")
    mkdir(opts.outputDir);
end

% Build summary table.
summaryTable = BuildCRLBSummaryTable(crlbLongTable, selectedMetabs, opts);

% Plot histograms.
for mIdx = 1:numel(selectedMetabs)

    metab = selectedMetabs(mIdx);

    vals = crlbLongTable.CRLB(string(crlbLongTable.metabolite) == metab);
    vals = vals(isfinite(vals));

    if isempty(vals)
        continue;
    end

    fig = figure( ...
        "Color", "w", ...
        "Name", "CRLB histogram - " + metab);

    histogram(vals, "BinWidth", opts.binWidth);

    hold on;

    xline(opts.crlbLimitLine, "r--", ...
        "CRLB = " + string(opts.crlbLimitLine), ...
        "LineWidth", 2, ...
        "LabelOrientation", "horizontal");

    xlabel("CRLB (%)");
    ylabel("Count");

    title("CRLB distribution: " + metab, ...
        "Interpreter", "none", ...
        "FontWeight", "bold");

    nTotal = numel(vals);
    pctUnderLimit = 100 * sum(vals < opts.crlbLimitLine) / nTotal;

    subtitle(sprintf("N = %d, %.1f%% with CRLB < %.0f", ...
        nTotal, pctUnderLimit, opts.crlbLimitLine));

    set(gca, ...
        "FontSize", 14, ...
        "FontWeight", "bold", ...
        "TickLabelInterpreter", "none", ...
        "XColor", [0 0 0], ...
        "YColor", [0 0 0]);

    grid on;

    if opts.saveFigures
        fileName = "CRLB_hist_" + matlab.lang.makeValidName(char(metab)) + ".png";
        exportgraphics(fig, fullfile(opts.outputDir, fileName), ...
            "Resolution", 300, ...
            "BackgroundColor", "white");
    end
end

% Export tables.
if opts.saveFigures
    writetable(crlbLongTable, fullfile(opts.outputDir, "CRLB_Long_Table.csv"));
    writetable(summaryTable, fullfile(opts.outputDir, "CRLB_Summary_Table.csv"));
end

crlbOutputs = struct();
crlbOutputs.longTable = crlbLongTable;
crlbOutputs.summaryTable = summaryTable;
crlbOutputs.selectedMetabolites = selectedMetabs;
crlbOutputs.options = opts;

end

% -------------------------------------------------------------------------
function opts = ParseCRLBHistogramOptions(cfg)

opts = struct();

opts.metabolites = GetOption(cfg, "metabolites", "all");

opts.useSumPreferredMetabolites = logical(GetOption(cfg, "useSumPreferredMetabolites", true));
opts.sumMetabolites = string(GetOption(cfg, "sumMetabolites", ...
    ["GPC+PCh", "NAA+NAAG", "Cr+PCr", "Glu+Gln"]));

opts.binWidth = double(GetOption(cfg, "binWidth", 5));
opts.crlbLimitLine = double(GetOption(cfg, "crlbLimitLine", 100));

opts.saveFigures = logical(GetOption(cfg, "saveFigures", true));
opts.outputDir = char(GetOption(cfg, "outputDir", fullfile(pwd, "CRLB_Histograms")));

end

% -------------------------------------------------------------------------
function value = GetOption(s, fieldName, defaultValue)

if isstruct(s) && isfield(s, fieldName)
    value = s.(fieldName);
else
    value = defaultValue;
end

end

% -------------------------------------------------------------------------
function patientTable = ExtractPatientCRLBLongTable(patientResult, patientID)

% This function tries several likely field names.
% The goal is to find a table containing raw part-level CRLB values.

patientID = string(patientID);

possibleFields = [ ...
    "crlbLongTable", ...
    "partCRLBTable", ...
    "partCrlbTable", ...
    "crlbTable", ...
    "CRLBTable" ...
];

for fIdx = 1:numel(possibleFields)

    fieldName = possibleFields(fIdx);

    if isfield(patientResult, fieldName)

        T = patientResult.(fieldName);

        if istable(T)
            patientTable = NormalizeCRLBTable(T, patientID);
            return;
        end
    end
end

error("Could not find raw part-level CRLB table for patient %s.", patientID);

end

% -------------------------------------------------------------------------
function Tlong = NormalizeCRLBTable(T, patientID)

% Expected output:
%   patientID | part | metabolite | CRLB

varNames = string(T.Properties.VariableNames);

% Case 1:
% Already long format.
if any(strcmpi(varNames, "metabolite")) && any(strcmpi(varNames, "CRLB"))

    metabCol = FindVariableIgnoreCase(T, "metabolite");
    crlbCol = FindVariableIgnoreCase(T, "CRLB");

    Tlong = table();
    Tlong.patientID = repmat(string(patientID), height(T), 1);

    if any(strcmpi(varNames, "patientID"))
        patientCol = FindVariableIgnoreCase(T, "patientID");
        Tlong.patientID = string(T{:, patientCol});
    end

    if any(strcmpi(varNames, "part"))
        partCol = FindVariableIgnoreCase(T, "part");
        Tlong.part = double(T{:, partCol});
    else
        Tlong.part = nan(height(T), 1);
    end

    Tlong.metabolite = string(T{:, metabCol});
    Tlong.CRLB = double(T{:, crlbCol});

    return;
end

% Case 2:
% Wide format: one row per part, metabolite names are columns.
% Example columns:
%   part | NAA | NAAG | NAA_NAAG | Cr | PCr | Cr_PCr | ...
bookkeeping = ["part", "patient", "patientID", "division", "file", "filename", "printFile", "coordFile"];
metabVars = strings(0, 1);

for k = 1:numel(varNames)
    if ~any(strcmpi(varNames(k), bookkeeping))
        metabVars(end + 1, 1) = varNames(k); %#ok<AGROW>
    end
end

if isempty(metabVars)
    error("Could not identify metabolite CRLB columns.");
end

if any(strcmpi(varNames, "part"))
    partCol = FindVariableIgnoreCase(T, "part");
    parts = double(T{:, partCol});
else
    parts = (1:height(T)).';
end

Tlong = table();

for mIdx = 1:numel(metabVars)

    metabVar = metabVars(mIdx);
    vals = double(T{:, metabVar});

    tmp = table();
    tmp.patientID = repmat(string(patientID), height(T), 1);
    tmp.part = parts(:);

    % Convert MATLAB-safe names like NAA_NAAG back to NAA+NAAG when possible.
    tmp.metabolite = repmat(PrettyMetaboliteName(metabVar), height(T), 1);
    tmp.CRLB = vals(:);

    Tlong = [Tlong; tmp]; %#ok<AGROW>
end

end

% -------------------------------------------------------------------------
function idx = FindVariableIgnoreCase(T, name)

varNames = string(T.Properties.VariableNames);
idx = find(strcmpi(varNames, string(name)), 1, "first");

if isempty(idx)
    error("Could not find variable %s.", string(name));
end

end

% -------------------------------------------------------------------------
function nameOut = PrettyMetaboliteName(nameIn)

nameOut = string(nameIn);

mapFrom = ["GPC_PCh", "NAA_NAAG", "Cr_PCr", "Glu_Gln"];
mapTo   = ["GPC+PCh", "NAA+NAAG", "Cr+PCr", "Glu+Gln"];

for k = 1:numel(mapFrom)
    if strcmp(nameOut, mapFrom(k))
        nameOut = mapTo(k);
        return;
    end
end

end

% -------------------------------------------------------------------------
function selected = ApplySumPreferredFiltering(metabsIn, sumMetabolites)

selected = string(metabsIn(:));
sumMetabolites = string(sumMetabolites(:));

for sIdx = 1:numel(sumMetabolites)

    sumName = sumMetabolites(sIdx);

    % Only remove components if the summed metabolite exists.
    if ~any(NamesMatchVector(selected, sumName))
        continue;
    end

    components = string(split(sumName, "+"));

    remove = false(numel(selected), 1);

    for cIdx = 1:numel(components)
        remove = remove | NamesMatchVector(selected, components(cIdx));
    end

    selected = selected(~remove);
end

end

% -------------------------------------------------------------------------
function tf = NamesMatchVector(names, target)

names = string(names(:));
target = string(target);

tf = strcmp(names, target);

targetValid = string(matlab.lang.makeValidName(char(target)));

for k = 1:numel(names)
    nameValid = string(matlab.lang.makeValidName(char(names(k))));
    tf(k) = tf(k) || strcmp(nameValid, targetValid);
end

end

% -------------------------------------------------------------------------
function summaryTable = BuildCRLBSummaryTable(crlbLongTable, selectedMetabs, opts)

summaryTable = table();

for mIdx = 1:numel(selectedMetabs)

    metab = selectedMetabs(mIdx);
    vals = crlbLongTable.CRLB(string(crlbLongTable.metabolite) == metab);
    vals = vals(isfinite(vals));

    row = table();

    row.metabolite = metab;
    row.nValues = numel(vals);

    if isempty(vals)
        row.meanCRLB = NaN;
        row.medianCRLB = NaN;
        row.minCRLB = NaN;
        row.maxCRLB = NaN;
        row.percentUnder20 = NaN;
        row.percentUnder30 = NaN;
        row.percentUnder100 = NaN;
        row.percentOverOrEqual100 = NaN;
    else
        row.meanCRLB = mean(vals);
        row.medianCRLB = median(vals);
        row.minCRLB = min(vals);
        row.maxCRLB = max(vals);
        row.percentUnder20 = 100 * sum(vals < 20) / numel(vals);
        row.percentUnder30 = 100 * sum(vals < 30) / numel(vals);
        row.percentUnder100 = 100 * sum(vals < opts.crlbLimitLine) / numel(vals);
        row.percentOverOrEqual100 = 100 * sum(vals >= opts.crlbLimitLine) / numel(vals);
    end

    summaryTable = [summaryTable; row]; %#ok<AGROW>
end

summaryTable = sortrows(summaryTable, "medianCRLB", "ascend");

end