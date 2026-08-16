function plotOutputs = PlotTemporalCorrelationDiagnostics(covOutputs, deGraafOutputs, plotCfg)
% PlotTemporalCorrelationDiagnostics
%
% Organized plotting and export pipeline for temporal metabolite-correlation analysis.
%
% This file wraps the exploratory plotting blocks into one callable function.
%
% Main features:
%   1. Empirical correlation matrix plots
%   2. LCModel / De Graaf amplitude-correlation matrix plots
%   3. Difference matrix: abs(LCModel) - abs(empirical)
%   4. CRLB reliability table and matrix overlay
%   5. Ranked empirical-vs-DeGraaf table with relative values and CRLB pass/fail
%   6. Per-patient z-scored time-series plots
%   7. Per-patient scatter plots with patient-specific CRLB status
%   8. Sum-metabolite scatter plots
%   9. Multi-patient pooled within-patient z-scored scatter plots
%
% Required inputs:
%   covOutputs     = MetabCovarianceByPatient(cfg)
%   deGraafOutputs = DeGraafAmplitudeCorrelationByPatient(cfg)
%
% Example call from RunExperiment.m:
%
%   plotCfg = struct();
%   plotCfg.patientIDs = "all";
%   plotCfg.matrixMetabs = "all";
%   plotCfg.scatterPairs = ["Glu", "GABA"; "Glu+Gln", "GABA"];
%   plotCfg.sumMetabs = ["GPC+PCh", "NAA+NAAG", "Cr+PCr", "Glu+Gln"];
%   plotCfg.saveFigures = true;
%   plotCfg.outputDir = fullfile(pwd, "TemporalCorrelationPlots");
%
%   plotOutputs = PlotTemporalCorrelationDiagnostics(covOutputs, deGraafOutputs, plotCfg);
%
% Notes:
%   - Patient-specific scatter plots use patient-specific CRLB reliability.
%   - Group matrix overlay and ranked table use group-level CRLB reliability.
%   - CRLB pass/fail is always calculated from fractionCRLBUnder100 at runtime.
%   - Pair passes CRLB only if BOTH metabolites pass the rule.

%% Defaults and setup

if nargin < 3
    plotCfg = struct();
end

plotCfg = ApplyPlotDefaults(plotCfg);

if plotCfg.saveFigures && ~isfolder(plotCfg.outputDir)
    mkdir(plotCfg.outputDir);
end

patientIDsAll = string(fieldnames(covOutputs.patientResultsByID));
patientIDs = ResolvePatientIDs(plotCfg.patientIDs, patientIDsAll);

fprintf('\nTemporal-correlation plotting pipeline\n');
fprintf('Selected patients: %d\n', numel(patientIDs));

plotOutputs = struct();
plotOutputs.selectedPatientIDs = patientIDs;
plotOutputs.settings = plotCfg;

%% Build group-level CRLB reliability table

groupCRLBQualityTable = BuildGroupCRLBQualityTable( ...
    covOutputs, ...
    plotCfg.crlbThreshold, ...
    plotCfg.requiredGoodFraction);

plotOutputs.groupCRLBQualityTable = groupCRLBQualityTable;

if plotCfg.exportCRLBQualityTable
    if ~isfolder(plotCfg.outputDir)
        mkdir(plotCfg.outputDir);
    end

    outputFile = fullfile(plotCfg.outputDir, "CRLB_reliability_table.csv");
    writetable(groupCRLBQualityTable, outputFile);
    fprintf('Saved CRLB reliability table to:\n%s\n', outputFile);
end

%% Matrix plots

if plotCfg.doMatrixPlots

    if isfield(covOutputs.group, 'meanCorrTable')
        PlotMatrixTable( ...
            covOutputs.group.meanCorrTable, ...
            "Group mean empirical correlation", ...
            [-1 1], ...
            plotCfg.matrixMetabs);

        MaybeSaveFigure(gcf, plotCfg, "group_mean_empirical_correlation");
    end

    if isfield(covOutputs.group, 'meanAbsCorrTable')
        PlotMatrixTable( ...
            covOutputs.group.meanAbsCorrTable, ...
            "Group mean absolute empirical correlation", ...
            [0 1], ...
            plotCfg.matrixMetabs);

        MaybeSaveFigure(gcf, plotCfg, "group_mean_abs_empirical_correlation");
    end

    if isfield(deGraafOutputs.group, 'meanAmplitudeCorrTable')
        PlotMatrixTable( ...
            deGraafOutputs.group.meanAmplitudeCorrTable, ...
            "Group mean LCModel / De Graaf amplitude correlation", ...
            [-1 1], ...
            plotCfg.matrixMetabs);

        MaybeSaveFigure(gcf, plotCfg, "group_mean_lcmodel_amplitude_correlation");
    end

    if isfield(deGraafOutputs.group, 'meanAbsAmplitudeCorrTable')
        PlotMatrixTable( ...
            deGraafOutputs.group.meanAbsAmplitudeCorrTable, ...
            "Group mean absolute LCModel / De Graaf amplitude correlation", ...
            [0 1], ...
            plotCfg.matrixMetabs);

        MaybeSaveFigure(gcf, plotCfg, "group_mean_abs_lcmodel_amplitude_correlation");
    end
end

%% Difference matrix

if plotCfg.doDifferenceMatrix

    T_emp_abs = GetRequiredGroupTable(covOutputs.group, 'meanAbsCorrTable');
    T_lcm_abs = GetRequiredGroupTable(deGraafOutputs.group, 'meanAbsAmplitudeCorrTable');

    PlotDifferenceMatrix( ...
        T_emp_abs, ...
        T_lcm_abs, ...
        plotCfg.matrixMetabs);

    MaybeSaveFigure(gcf, plotCfg, "difference_abs_lcmodel_minus_abs_empirical");
end

%% Empirical matrix with CRLB overlay

if plotCfg.doCRLBOverlayMatrix

    T_emp_abs = GetRequiredGroupTable(covOutputs.group, 'meanAbsCorrTable');

    PlotEmpiricalMatrixWithCRLBOverlay( ...
        T_emp_abs, ...
        groupCRLBQualityTable, ...
        plotCfg.requiredGoodFraction, ...
        plotCfg.matrixMetabs);

    MaybeSaveFigure(gcf, plotCfg, "mean_abs_empirical_with_crlb_overlay");
end

%% Ranked empirical-vs-LCModel table

if plotCfg.doRankedTable

    T_emp_abs = GetRequiredGroupTable(covOutputs.group, 'meanAbsCorrTable');
    T_lcm_abs = GetRequiredGroupTable(deGraafOutputs.group, 'meanAbsAmplitudeCorrTable');

    rankedTable = BuildRankedEmpiricalDeGraafTable( ...
        T_emp_abs, ...
        T_lcm_abs, ...
        groupCRLBQualityTable, ...
        plotCfg.requiredGoodFraction, ...
        plotCfg.matrixMetabs);

    plotOutputs.rankedEmpiricalDeGraafTable = rankedTable;

    if plotCfg.exportRankedTable
        if ~isfolder(plotCfg.outputDir)
            mkdir(plotCfg.outputDir);
        end

        outputFile = fullfile(plotCfg.outputDir, ...
            "Ranked_Empirical_vs_DeGraaf_with_CRLB.csv");

        writetable(rankedTable, outputFile);
        fprintf('Saved ranked table to:\n%s\n', outputFile);
    end
end

%% Per-patient z-scored time-series plots

if plotCfg.doTimeSeriesPlots

    for pIdx = 1:numel(patientIDs)

        patientID = patientIDs(pIdx);
        safePatientID = matlab.lang.makeValidName(char(patientID));
        T = covOutputs.patientResultsByID.(safePatientID).partTable;

        PlotTimeSeriesGroupsForPatient( ...
            T, ...
            patientID, ...
            plotCfg.timeSeriesGroups);

        MaybeSaveFigure(gcf, plotCfg, "time_series_groups_" + patientID);
    end
end

%% Per-patient scatter plots

if plotCfg.doPatientScatterPlots || plotCfg.doSumScatterPlots

    for pIdx = 1:numel(patientIDs)

        patientID = patientIDs(pIdx);
        safePatientID = matlab.lang.makeValidName(char(patientID));

        T = covOutputs.patientResultsByID.(safePatientID).partTable;
        coordTable = covOutputs.patientResultsByID.(safePatientID).coordTable;

        patientCRLBQualityTable = BuildPatientCRLBQualityTable( ...
            covOutputs, ...
            patientID, ...
            coordTable, ...
            T, ...
            plotCfg.crlbThreshold, ...
            plotCfg.requiredGoodFraction);

        if plotCfg.doPatientScatterPlots
            PlotScatterPairsForPatient( ...
                T, ...
                patientID, ...
                plotCfg.scatterPairs, ...
                patientCRLBQualityTable, ...
                plotCfg.requiredGoodFraction);

            MaybeSaveFigure(gcf, plotCfg, "scatter_pairs_" + patientID);
        end

        if plotCfg.doSumScatterPlots
            PlotSumScatterForPatient( ...
                T, ...
                patientID, ...
                plotCfg.sumMetabs, ...
                patientCRLBQualityTable, ...
                plotCfg.requiredGoodFraction);

            MaybeSaveFigure(gcf, plotCfg, "sum_scatter_" + patientID);
        end
    end
end

%% Multi-patient pooled within-patient z-scored scatter

if plotCfg.doPooledZScatter

    for pairIdx = 1:size(plotCfg.pooledZPairs, 1)

        xName = plotCfg.pooledZPairs(pairIdx, 1);
        yName = plotCfg.pooledZPairs(pairIdx, 2);

        PlotPooledZScatter( ...
            covOutputs, ...
            patientIDs, ...
            xName, ...
            yName);

        MaybeSaveFigure(gcf, plotCfg, "pooled_z_scatter_" + xName + "_vs_" + yName);
    end
end

fprintf('Finished temporal-correlation plotting pipeline.\n');

end

%% ========================================================================
% Defaults and options
% ========================================================================

function plotCfg = ApplyPlotDefaults(plotCfg)

if nargin < 1 || isempty(plotCfg)
    plotCfg = struct();
end

plotCfg = SetDefault(plotCfg, 'patientIDs', "all");

plotCfg = SetDefault(plotCfg, 'doMatrixPlots', true);
plotCfg = SetDefault(plotCfg, 'doDifferenceMatrix', true);
plotCfg = SetDefault(plotCfg, 'doCRLBOverlayMatrix', true);
plotCfg = SetDefault(plotCfg, 'matrixMetabs', "all");

plotCfg = SetDefault(plotCfg, 'crlbThreshold', 100);
plotCfg = SetDefault(plotCfg, 'requiredGoodFraction', 0.90);

plotCfg = SetDefault(plotCfg, 'doTimeSeriesPlots', true);
plotCfg = SetDefault(plotCfg, 'timeSeriesGroups', { ...
    ["GPC", "PCh", "GPC+PCh"], ...
    ["Glu", "Gln", "Glu+Gln"], ...
    ["NAA", "NAAG", "NAA+NAAG"], ...
    ["Cr", "PCr", "Cr+PCr"]});

plotCfg = SetDefault(plotCfg, 'doPatientScatterPlots', true);
plotCfg = SetDefault(plotCfg, 'scatterPairs', ["Glu", "GABA"]);

plotCfg = SetDefault(plotCfg, 'doSumScatterPlots', true);
plotCfg = SetDefault(plotCfg, 'sumMetabs', ...
    ["GPC+PCh", "NAA+NAAG", "Cr+PCr", "Glu+Gln"]);

plotCfg = SetDefault(plotCfg, 'doPooledZScatter', true);
plotCfg = SetDefault(plotCfg, 'pooledZPairs', ["Glu", "GABA"]);

plotCfg = SetDefault(plotCfg, 'doRankedTable', true);
plotCfg = SetDefault(plotCfg, 'exportRankedTable', true);
plotCfg = SetDefault(plotCfg, 'exportCRLBQualityTable', true);

plotCfg = SetDefault(plotCfg, 'saveFigures', false);
plotCfg = SetDefault(plotCfg, 'outputDir', fullfile(pwd, "TemporalCorrelationPlots"));

end

function s = SetDefault(s, fieldName, defaultValue)

if ~isfield(s, fieldName) || isempty(s.(fieldName))
    s.(fieldName) = defaultValue;
end

end

function patientIDs = ResolvePatientIDs(patientOption, patientIDsAll)

if ischar(patientOption) || isstring(patientOption)

    patientOption = string(patientOption);

    if isscalar(patientOption) && patientOption == "all"
        patientIDs = patientIDsAll;
    else
        patientIDs = patientOption(:);
    end

elseif isnumeric(patientOption)

    patientIDs = patientIDsAll(patientOption);

else

    error('Unsupported plotCfg.patientIDs format. Use "all", string array, or numeric indices.');
end

missing = patientIDs(~ismember(patientIDs, patientIDsAll));

if ~isempty(missing)
    error('Requested patient IDs were not found: %s', strjoin(missing, ", "));
end

end

%% ========================================================================
% Matrix plotting
% ========================================================================

function PlotMatrixTable(T, plotTitle, climits, selectedMetabs)

[Tsub, labels] = SubsetSquareTable(T, selectedMetabs);
M = table2array(Tsub);

figure;
imagesc(M);
axis image;
colorbar;

caxis(climits);
colormap(gray);

title(plotTitle, 'Interpreter', 'none');

xticks(1:numel(labels));
xticklabels(labels);
xtickangle(45);

yticks(1:numel(labels));
yticklabels(labels);

set(gca, 'TickLabelInterpreter', 'none');
set(gca, 'FontSize', 12);

end

function PlotDifferenceMatrix(T_emp, T_lcm, selectedMetabs)

[T_emp_sub, labelsEmp] = SubsetSquareTable(T_emp, selectedMetabs);
[T_lcm_sub, labelsLCM] = SubsetSquareTable(T_lcm, selectedMetabs);

commonLabels = labelsEmp(ismember(labelsEmp, labelsLCM));

[~, idxEmp] = ismember(commonLabels, labelsEmp);
[~, idxLCM] = ismember(commonLabels, labelsLCM);

M_emp = table2array(T_emp_sub(idxEmp, idxEmp));
M_lcm = table2array(T_lcm_sub(idxLCM, idxLCM));

D = M_lcm - M_emp;

figure;
imagesc(D);
axis image;
colorbar;

lim = max(abs(D(:)), [], 'omitnan');

if isempty(lim) || isnan(lim) || lim == 0
    lim = 1;
end

caxis([-lim lim]);
colormap(gray);

title("Difference: abs LCModel / De Graaf - abs empirical", ...
    'Interpreter', 'none');

xticks(1:numel(commonLabels));
xticklabels(commonLabels);
xtickangle(45);

yticks(1:numel(commonLabels));
yticklabels(commonLabels);

set(gca, 'TickLabelInterpreter', 'none');
set(gca, 'FontSize', 12);

end

function PlotEmpiricalMatrixWithCRLBOverlay(Tcorr, crlbQualityTable, requiredGoodFraction, selectedMetabs)

[Tsub, labels] = SubsetSquareTable(Tcorr, selectedMetabs);
M = table2array(Tsub);

crlbQualityTable.metabolite = string(crlbQualityTable.metabolite);

[tf, idx] = ismember(labels, crlbQualityTable.metabolite);

if any(~tf)
    missing = labels(~tf);
    error('Some matrix labels were not found in CRLB table: %s', strjoin(missing, ", "));
end

frac = crlbQualityTable.fractionCRLBUnder100(idx);
badMetabMask = frac < requiredGoodFraction;
badCellMask = badMetabMask(:) | badMetabMask(:).';

figure;
imagesc(M);
axis image;
colorbar;

caxis([0 1]);
colormap(gray);

title("Mean absolute empirical correlation with CRLB reliability warning", ...
    'Interpreter', 'none');

xticks(1:numel(labels));
xticklabels(labels);
xtickangle(45);

yticks(1:numel(labels));
yticklabels(labels);

set(gca, 'TickLabelInterpreter', 'none');
set(gca, 'FontSize', 12);

hold on;

[rowBad, colBad] = find(badCellMask);

plot(colBad, rowBad, 'wo', ...
    'MarkerSize', 8, ...
    'LineWidth', 1.2, ...
    'HandleVisibility', 'off');

plot(colBad, rowBad, 'kx', ...
    'MarkerSize', 7, ...
    'LineWidth', 1.4, ...
    'DisplayName', 'CRLB reliability fail');

hold off;

legend('Location', 'eastoutside');

end

function [Tsub, labels] = SubsetSquareTable(T, selectedMetabs)

labelsAll = string(T.Properties.RowNames);

if ischar(selectedMetabs) || isstring(selectedMetabs)

    selectedMetabs = string(selectedMetabs);

    if isscalar(selectedMetabs) && selectedMetabs == "all"
        labels = labelsAll;
    else
        labels = selectedMetabs(:);
    end

else

    labels = string(selectedMetabs(:));
end

missing = labels(~ismember(labels, labelsAll));

if ~isempty(missing)
    error('Requested metabolites were not found in matrix table: %s', strjoin(missing, ", "));
end

[~, idx] = ismember(labels, labelsAll);
Tsub = T(idx, idx);

end

function T = GetRequiredGroupTable(groupStruct, fieldName)

if ~isfield(groupStruct, fieldName)
    error('Required group table was not found: %s', fieldName);
end

T = groupStruct.(fieldName);

end

%% ========================================================================
% CRLB reliability
% ========================================================================

function crlbQualityTable = BuildGroupCRLBQualityTable(covOutputs, crlbThreshold, requiredGoodFraction)

metabList = string(covOutputs.metabList(:));
nMetabs = numel(metabList);

goodCRLBCount = zeros(nMetabs, 1);
totalInstanceCount = zeros(nMetabs, 1);

patientIDs = string(fieldnames(covOutputs.patientResultsByID));

for pIdx = 1:numel(patientIDs)

    patientID = patientIDs(pIdx);
    safePatientID = matlab.lang.makeValidName(char(patientID));

    Tparts = covOutputs.patientResultsByID.(safePatientID).partTable;
    Tcoord = covOutputs.patientResultsByID.(safePatientID).coordTable;

    patientTable = BuildPatientCRLBQualityTable( ...
        covOutputs, ...
        patientID, ...
        Tcoord, ...
        Tparts, ...
        crlbThreshold, ...
        requiredGoodFraction);

    [tf, idx] = ismember(metabList, patientTable.metabolite);

    goodCRLBCount(tf) = goodCRLBCount(tf) + patientTable.nCRLBUnder100(idx(tf));
    totalInstanceCount(tf) = totalInstanceCount(tf) + patientTable.nInstances(idx(tf));
end

fractionCRLBUnder100 = goodCRLBCount ./ totalInstanceCount;
fails90PercentRule = fractionCRLBUnder100 < requiredGoodFraction;

crlbQualityTable = table( ...
    metabList, ...
    goodCRLBCount, ...
    totalInstanceCount, ...
    fractionCRLBUnder100, ...
    fails90PercentRule, ...
    'VariableNames', { ...
    'metabolite', ...
    'nCRLBUnder100', ...
    'nInstances', ...
    'fractionCRLBUnder100', ...
    'fails90PercentRule'});

end

function patientCRLBQualityTable = BuildPatientCRLBQualityTable(covOutputs, patientID, coordTable, partTable, crlbThreshold, requiredGoodFraction)

metabList = string(covOutputs.metabList(:));
nMetabs = numel(metabList);

if isfield(covOutputs, 'settings') && isfield(covOutputs.settings, 'division')
    divisionUsed = covOutputs.settings.division;
else
    divisionUsed = 1;
end

partsUsed = GetPartVector(partTable);

coordTable = EnsureCoordTableStandardColumns(coordTable);
crlbCol = FindCRLBColumn(coordTable);

nRows = height(coordTable);
parsedDivision = nan(nRows, 1);
parsedPart = nan(nRows, 1);

for r = 1:nRows

    [~, baseName, ext] = fileparts(coordTable.filename(r));
    curFile = string(baseName) + string(ext);

    tok = regexp(curFile, ...
        '.*Division_(\d+)_(?:part_)?(\d+)\.basis\.coord$', ...
        'tokens', 'once');

    if isempty(tok)
        continue;
    end

    parsedDivision(r) = str2double(tok{1});
    parsedPart(r) = str2double(tok{2});
end

coordTable.division = parsedDivision;
coordTable.part = parsedPart;
coordTable = coordTable(coordTable.division == divisionUsed, :);

goodCRLBCount = zeros(nMetabs, 1);
totalInstanceCount = zeros(nMetabs, 1);
fractionCRLBUnder100 = nan(nMetabs, 1);
fails90PercentRule = false(nMetabs, 1);

for m = 1:nMetabs

    metabName = metabList(m);
    crlbVals = nan(numel(partsUsed), 1);

    for partIdx = 1:numel(partsUsed)

        curPart = partsUsed(partIdx);
        idx = coordTable.name == metabName & coordTable.part == curPart;

        if sum(idx) >= 1
            tmp = coordTable.(char(crlbCol))(idx);
            crlbVals(partIdx) = double(tmp(1));
        end
    end

    totalInstanceCount(m) = numel(crlbVals);
    goodCRLBCount(m) = sum(crlbVals < crlbThreshold, 'omitnan');

    fractionCRLBUnder100(m) = goodCRLBCount(m) / totalInstanceCount(m);
    fails90PercentRule(m) = fractionCRLBUnder100(m) < requiredGoodFraction;
end

patientCRLBQualityTable = table( ...
    metabList, ...
    goodCRLBCount, ...
    totalInstanceCount, ...
    fractionCRLBUnder100, ...
    fails90PercentRule, ...
    'VariableNames', { ...
    'metabolite', ...
    'nCRLBUnder100', ...
    'nInstances', ...
    'fractionCRLBUnder100', ...
    'fails90PercentRule'});

end

function coordTable = EnsureCoordTableStandardColumns(coordTable)

varNames = string(coordTable.Properties.VariableNames);

nameCandidates = ["name", "metab", "metabolite", "metabName", "Metabolite"];
nameCol = "";

for c = nameCandidates
    if ismember(c, varNames)
        nameCol = c;
        break;
    end
end

if strlength(nameCol) == 0
    error('Could not find a metabolite-name column in coordTable. Available columns are: %s', strjoin(varNames, ", "));
end

filenameCandidates = ["filename", "fileName", "FileName", "coordFile", "file"];
filenameCol = "";

for c = filenameCandidates
    if ismember(c, varNames)
        filenameCol = c;
        break;
    end
end

if strlength(filenameCol) == 0
    error('Could not find a filename column in coordTable. Available columns are: %s', strjoin(varNames, ", "));
end

coordTable.name = string(coordTable.(char(nameCol)));
coordTable.filename = string(coordTable.(char(filenameCol)));

end

function crlbCol = FindCRLBColumn(T)

varNames = string(T.Properties.VariableNames);

crlbCandidates = ["CRLB", "crlb", "SD", "sd", ...
    "percentSD", "PercentSD", "pctSD", "pctCrLB", ...
    "crlbPercent", "CRLBPercent"];

crlbCol = "";

for c = crlbCandidates
    if ismember(c, varNames)
        crlbCol = c;
        return;
    end
end

error('Could not find a CRLB / %%SD column. Available columns are: %s', ...
    strjoin(varNames, ", "));

end

function parts = GetPartVector(partTable)

if ismember("part", string(partTable.Properties.VariableNames))
    parts = partTable.part;
else
    parts = (1:height(partTable)).';
end

end

%% ========================================================================
% Ranked table
% ========================================================================

function rankedTable = BuildRankedEmpiricalDeGraafTable(T_emp, T_lcm, crlbQualityTable, requiredGoodFraction, selectedMetabs)

[T_emp_sub, labelsEmp] = SubsetSquareTable(T_emp, selectedMetabs);
[T_lcm_sub, labelsLCM] = SubsetSquareTable(T_lcm, selectedMetabs);

commonLabels = labelsEmp(ismember(labelsEmp, labelsLCM));

[~, idxEmp] = ismember(commonLabels, labelsEmp);
[~, idxLCM] = ismember(commonLabels, labelsLCM);

M_emp = table2array(T_emp_sub(idxEmp, idxEmp));
M_lcm = table2array(T_lcm_sub(idxLCM, idxLCM));

diff_LCM_minus_empirical = M_lcm - M_emp;
absDifference = abs(diff_LCM_minus_empirical);

upperMask = triu(true(numel(commonLabels)), 1);

[rowIdx, colIdx] = find(upperMask);

metaboliteA = commonLabels(rowIdx);
metaboliteB = commonLabels(colIdx);

absEmpiricalCorr = M_emp(upperMask);
absLCModelCorr = M_lcm(upperMask);

diffVals = diff_LCM_minus_empirical(upperMask);
absDiffVals = absDifference(upperMask);

rankedTable = table( ...
    metaboliteA, ...
    metaboliteB, ...
    absEmpiricalCorr, ...
    absLCModelCorr, ...
    diffVals, ...
    absDiffVals, ...
    'VariableNames', { ...
    'metaboliteA', ...
    'metaboliteB', ...
    'absEmpiricalCorr', ...
    'absLCModelCorr', ...
    'diff_LCM_minus_empirical', ...
    'absDifference'});

rankedTable = rankedTable(~isnan(rankedTable.diff_LCM_minus_empirical), :);

emp = rankedTable.absEmpiricalCorr;
lcm = rankedTable.absLCModelCorr;

empirical_over_LCModel = nan(height(rankedTable), 1);
LCModel_over_empirical = nan(height(rankedTable), 1);
relativeDiff_LCM_minus_empirical_percent = nan(height(rankedTable), 1);

validLCM = ~isnan(lcm) & lcm ~= 0;
validEmp = ~isnan(emp) & emp ~= 0;

empirical_over_LCModel(validLCM) = emp(validLCM) ./ lcm(validLCM);
LCModel_over_empirical(validEmp) = lcm(validEmp) ./ emp(validEmp);

relativeDiff_LCM_minus_empirical_percent(validLCM) = ...
    100 * (lcm(validLCM) - emp(validLCM)) ./ lcm(validLCM);

rankedTable.empirical_over_LCModel = empirical_over_LCModel;
rankedTable.LCModel_over_empirical = LCModel_over_empirical;
rankedTable.relativeDiff_LCM_minus_empirical_percent = relativeDiff_LCM_minus_empirical_percent;

crlbQualityTable.metabolite = string(crlbQualityTable.metabolite);

[tfA, idxA] = ismember(rankedTable.metaboliteA, crlbQualityTable.metabolite);
[tfB, idxB] = ismember(rankedTable.metaboliteB, crlbQualityTable.metabolite);

if any(~tfA)
    missingA = unique(rankedTable.metaboliteA(~tfA));
    error('Some metaboliteA names were not found in CRLB table: %s', strjoin(missingA, ", "));
end

if any(~tfB)
    missingB = unique(rankedTable.metaboliteB(~tfB));
    error('Some metaboliteB names were not found in CRLB table: %s', strjoin(missingB, ", "));
end

fracA = crlbQualityTable.fractionCRLBUnder100(idxA);
fracB = crlbQualityTable.fractionCRLBUnder100(idxB);

passA = fracA >= requiredGoodFraction;
passB = fracB >= requiredGoodFraction;

pairPassesCRLB = passA & passB;

CRLB_pair_status = strings(height(rankedTable), 1);
CRLB_pair_status(pairPassesCRLB) = "PASS";
CRLB_pair_status(~pairPassesCRLB) = "FAIL";

rankedTable.CRLB_pair_status = CRLB_pair_status;
rankedTable.metaboliteA_fractionCRLBUnder100 = fracA;
rankedTable.metaboliteB_fractionCRLBUnder100 = fracB;
rankedTable.metaboliteA_failsCRLB = ~passA;
rankedTable.metaboliteB_failsCRLB = ~passB;

rankedTable = sortrows(rankedTable, "absDifference", "descend");

end

%% ========================================================================
% Time-series plotting
% ========================================================================

function PlotTimeSeriesGroupsForPatient(T, patientID, timeSeriesGroups)

nGroups = numel(timeSeriesGroups);

nCols = ceil(sqrt(nGroups));
nRows = ceil(nGroups / nCols);

figure;
tiledlayout(nRows, nCols, ...
    'TileSpacing', 'compact', ...
    'Padding', 'compact');

for g = 1:nGroups

    metabs = string(timeSeriesGroups{g});
    nexttile;
    hold on;

    for m = 1:numel(metabs)

        metabName = metabs(m);
        colName = matlab.lang.makeValidName(char(metabName));

        if ~ismember(colName, string(T.Properties.VariableNames))
            warning('Patient %s: metabolite %s was not found in partTable.', patientID, metabName);
            continue;
        end

        y = T.(colName);
        parts = GetPartVector(T);
        z = ZScoreOmitNaN(y);

        plot(parts, z, '-o', 'DisplayName', metabName);
    end

    yline(0, '--');
    hold off;
    grid on;

    xlabel('Division_1 part');
    ylabel('Within-patient z-score');

    title(strjoin(metabs, ", "), 'Interpreter', 'none');

    set(gca, 'FontSize', 9);
end

legend('Location', 'eastoutside');

sgtitle("Z-scored metabolite time-series - " + patientID, ...
    'Interpreter', 'none');

end

%% ========================================================================
% Scatter plotting
% ========================================================================

function PlotScatterPairsForPatient(T, patientID, scatterPairs, patientCRLBQualityTable, requiredGoodFraction)

if isempty(scatterPairs)
    return;
end

nPairs = size(scatterPairs, 1);

nCols = ceil(sqrt(nPairs));
nRows = ceil(nPairs / nCols);

figure;
tiledlayout(nRows, nCols, ...
    'TileSpacing', 'compact', ...
    'Padding', 'compact');

for pairIdx = 1:nPairs

    xName = scatterPairs(pairIdx, 1);
    yName = scatterPairs(pairIdx, 2);

    nexttile;

    PlotSingleScatterInCurrentAxes( ...
        T, ...
        patientID, ...
        xName, ...
        yName, ...
        patientCRLBQualityTable, ...
        requiredGoodFraction);
end

sgtitle("Selected scatter plots with patient-specific CRLB status - " + patientID, ...
    'Interpreter', 'none');

end

function PlotSumScatterForPatient(T, patientID, sumMetabs, patientCRLBQualityTable, requiredGoodFraction)

sumMetabs = string(sumMetabs(:));
nSums = numel(sumMetabs);

if nSums < 2
    return;
end

nPairs = nchoosek(nSums, 2);

nCols = ceil(sqrt(nPairs));
nRows = ceil(nPairs / nCols);

figure;
tiledlayout(nRows, nCols, ...
    'TileSpacing', 'compact', ...
    'Padding', 'compact');

for a = 1:nSums-1

    for b = a+1:nSums

        xName = sumMetabs(a);
        yName = sumMetabs(b);

        nexttile;

        PlotSingleScatterInCurrentAxes( ...
            T, ...
            patientID, ...
            xName, ...
            yName, ...
            patientCRLBQualityTable, ...
            requiredGoodFraction);
    end
end

sgtitle("Sum-metabolite scatter plots with patient-specific CRLB status - " + patientID, ...
    'Interpreter', 'none');

end

function PlotSingleScatterInCurrentAxes(T, patientID, xName, yName, patientCRLBQualityTable, requiredGoodFraction)

xCol = matlab.lang.makeValidName(char(xName));
yCol = matlab.lang.makeValidName(char(yName));

if ~ismember(xCol, string(T.Properties.VariableNames)) || ...
   ~ismember(yCol, string(T.Properties.VariableNames))

    text(0.5, 0.5, ...
        sprintf('Missing data:\n%s vs %s', yName, xName), ...
        'HorizontalAlignment', 'center', ...
        'Interpreter', 'none');

    axis off;
    return;
end

x = T.(xCol);
y = T.(yCol);

valid = ~isnan(x) & ~isnan(y);
nValid = sum(valid);

if nValid >= 3
    rVal = corr(x(valid), y(valid));
else
    rVal = NaN;
end

[fracX, fracY, crlbStatus] = GetPairCRLBStatus( ...
    xName, ...
    yName, ...
    patientCRLBQualityTable, ...
    requiredGoodFraction);

scatter(x(valid), y(valid), 45, 'filled');
grid on;

xlabel(xName, 'Interpreter', 'none');
ylabel(yName, 'Interpreter', 'none');

title(sprintf('%s vs %s\nr = %.2f, n = %d | %s\n%s %.0f%%, %s %.0f%% CRLB<100', ...
    yName, xName, ...
    rVal, nValid, ...
    crlbStatus, ...
    xName, 100 * fracX, ...
    yName, 100 * fracY), ...
    'Interpreter', 'none');

set(gca, 'FontSize', 9);

end

function [fracX, fracY, crlbStatus] = GetPairCRLBStatus(xName, yName, crlbTable, requiredGoodFraction)

crlbTable.metabolite = string(crlbTable.metabolite);

idxX = find(crlbTable.metabolite == string(xName), 1);
idxY = find(crlbTable.metabolite == string(yName), 1);

if isempty(idxX)
    fracX = NaN;
else
    fracX = crlbTable.fractionCRLBUnder100(idxX);
end

if isempty(idxY)
    fracY = NaN;
else
    fracY = crlbTable.fractionCRLBUnder100(idxY);
end

passX = fracX >= requiredGoodFraction;
passY = fracY >= requiredGoodFraction;

if passX && passY
    crlbStatus = "CRLB PASS";
else
    crlbStatus = "CRLB FAIL";
end

end

%% ========================================================================
% Pooled within-patient z-scored scatter
% ========================================================================

function PlotPooledZScatter(covOutputs, patientIDs, xName, yName)

allX = [];
allY = [];
allPatient = strings(0, 1);

perPatientR = nan(numel(patientIDs), 1);
perPatientN = nan(numel(patientIDs), 1);

for pIdx = 1:numel(patientIDs)

    patientID = patientIDs(pIdx);
    safePatientID = matlab.lang.makeValidName(char(patientID));
    T = covOutputs.patientResultsByID.(safePatientID).partTable;

    xCol = matlab.lang.makeValidName(char(xName));
    yCol = matlab.lang.makeValidName(char(yName));

    if ~ismember(xCol, string(T.Properties.VariableNames)) || ...
       ~ismember(yCol, string(T.Properties.VariableNames))
        continue;
    end

    x = T.(xCol);
    y = T.(yCol);

    valid = ~isnan(x) & ~isnan(y);
    nValid = sum(valid);

    perPatientN(pIdx) = nValid;

    if nValid < 3
        continue;
    end

    xv = x(valid);
    yv = y(valid);

    perPatientR(pIdx) = corr(xv, yv);

    zx = ZScoreOmitNaN(xv);
    zy = ZScoreOmitNaN(yv);

    allX = [allX; zx]; %#ok<AGROW>
    allY = [allY; zy]; %#ok<AGROW>
    allPatient = [allPatient; repmat(patientID, numel(zx), 1)]; %#ok<AGROW>
end

validPool = ~isnan(allX) & ~isnan(allY);

if sum(validPool) >= 3
    pooledWithinPatientR = corr(allX(validPool), allY(validPool));
else
    pooledWithinPatientR = NaN;
end

meanPatientR = mean(perPatientR, 'omitnan');
medianPatientR = median(perPatientR, 'omitnan');

figure;
gscatter(allX, allY, allPatient);
grid on;

xlabel("Within-patient z-score: " + xName, 'Interpreter', 'none');
ylabel("Within-patient z-score: " + yName, 'Interpreter', 'none');

title(sprintf(['Multi-patient pooled scatter: %s vs %s\n' ...
    'pooled within-patient z r = %.2f | mean patient r = %.2f | median patient r = %.2f'], ...
    yName, xName, pooledWithinPatientR, meanPatientR, medianPatientR), ...
    'Interpreter', 'none');

legend('Location', 'eastoutside');
set(gca, 'FontSize', 12);

end

%% ========================================================================
% Utilities
% ========================================================================

function z = ZScoreOmitNaN(x)

mu = mean(x, 'omitnan');
sigma = std(x, 0, 'omitnan');

if isempty(sigma) || isnan(sigma) || sigma == 0
    z = nan(size(x));
else
    z = (x - mu) ./ sigma;
end

end

function MaybeSaveFigure(figHandle, plotCfg, fileBaseName)

if ~plotCfg.saveFigures
    return;
end

if ~isfolder(plotCfg.outputDir)
    mkdir(plotCfg.outputDir);
end

safeName = string(matlab.lang.makeValidName(char(fileBaseName)));
outputFile = fullfile(plotCfg.outputDir, safeName + ".png");

saveas(figHandle, outputFile);

fprintf('Saved figure: %s\n', outputFile);

end
