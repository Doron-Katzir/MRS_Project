%% 
clear;
clc;

cfg = ProjectConfig();

%% User input for this run

cfg.paths.rootDir = "C:\Users\User\Documents\Thesis_Lab";
cfg.paths.dataDir = fullfile(cfg.paths.rootDir, "Data");
cfg.paths.coordDir = fullfile(cfg.paths.rootDir, "LCMFit");

%% Splice data and save coord files
%% Splice only one file

% cfg.input.mode = "singleFile";
% cfg.input.singlePatientID = "P01";
% cfg.input.singleFileName = "meas_MID00020_FID54986_eja_svs_slaser_TE_80_1_PRE.dat";
% cfg.input.singleFile = fullfile(cfg.paths.dataDir, cfg.input.singleFileName);

%% Splice all files in the working directory

% cfg.input.mode = "directory";
% 
% cfg.input.directory = cfg.paths.dataDir;
% cfg.input.filePattern = "*.dat";
% cfg.input.recursive = false;

%% Call splice function
% processedPatients = Splice_data_multi_patient_input(cfg);

%% Load coord files
%% Load only one fitted patient folder

% cfg.load.mode = "singleSubfolder";
% cfg.load.selectedPatientID = "P11";
% cfg.load.coordDir = cfg.paths.coordDir;
% 

%% Load all subfolders

cfg.load.mode = "allSubfolders";
cfg.load.coordDir = cfg.paths.coordDir;

%% Check bias
% outputs = Load_subsets_multi_patient_input(cfg);


%% Check metab covariance

cfg.covariance.loadMode = "allSubfolders";
% Preserve raw zero values until the centralized analysis-filtering stage.
cfg.covariance.ignoreZeros = false;
covOutputs = MetabCovarianceByPatient(cfg);



%% Show average empirical metabolite correlation matrix as image

T = covOutputs.group.meanCorrTable;

M = table2array(T);
labels = string(T.Properties.RowNames);

figure;
imagesc(M);
axis image;
colorbar;

caxis([-1 1]);

colormap(gray);
title("Average empirical metabolite correlation matrix", ...
    'Interpreter', 'none');

xticks(1:numel(labels));
xticklabels(labels);
xtickangle(45);

yticks(1:numel(labels));
yticklabels(labels);

set(gca, 'TickLabelInterpreter', 'none');
set(gca, 'FontSize', 15);

%% De graaf Amplitude Correlation

%% All subfolders
cfg.degraaf.loadMode = "allSubfolders";
cfg.degraaf.division = 1;
cfg.degraaf.maskInvalidCRLB = false;
cfg.degraaf.invalidCRLBValue = 100;

%% Subset of subfolders
% cfg.degraaf.loadMode = "selectedSubfolders";
% cfg.degraaf.selectedPatientIDs = ["P01", "P03", "P11"];
% cfg.degraaf.division = 1;


%% Run

deGraafOutputs = DeGraafAmplitudeCorrelationByPatient(cfg);
% deGraafOutputs.patientResultsByID.P11.meanAmplitudeCovTable



%% Show De Graaf amplitude correlation matrix as image

T = deGraafOutputs.group.meanAmplitudeCorrTable;

M = table2array(T);
labels = string(T.Properties.RowNames);

figure;
imagesc(M);
axis image;
colorbar;

caxis([-1 1]);

colormap(gray);

xticks(1:numel(labels));
xticklabels(labels);
xtickangle(45);

yticks(1:numel(labels));
yticklabels(labels);

set(gca, 'TickLabelInterpreter', 'none');
set(gca, 'FontSize', 15);

title("Group mean De Graaf / LCModel amplitude correlation", ...
    'Interpreter', 'none');


%% Plotting configurations

plotCfg = struct();

% Patient selection
% Options:
%   "all"
%   ["P01", "P03", "P11"]
%   1:5
plotCfg.patientIDs = 1:10;

% Matrix plots
% Example subset:
% plotCfg.matrixMetabs = ["Glu", "Gln", "Glu+Gln", "GABA", "NAA+NAAG", "Cr+PCr"];
plotCfg.doMatrixPlots = false;
plotCfg.doDifferenceMatrix = false;
plotCfg.doCRLBOverlayMatrix = false;
plotCfg.matrixMetabs = "all";


% CRLB rule
plotCfg.crlbThreshold = 100;
plotCfg.requiredGoodFraction = 0.01;

% Time-series plots
plotCfg.doTimeSeriesPlots = true;
plotCfg.timeSeriesGroups = {
    ["NAA+NAAG","Glu_Gln"]};

% Per-patient scatter plots
plotCfg.doPatientScatterPlots = true;
plotCfg.scatterPairs = ["NAA+NAAG","Glu+Gln"];

% Sum-metabolite scatter plots
% plotCfg.doSumScatterPlots = true;
% plotCfg.sumMetabs = ["GPC+PCh", "NAA+NAAG", "Cr+PCr", "Glu+Gln"];

% Multi-patient pooled within-patient z-scored scatter
% plotCfg.doPooledZScatter = true;
% plotCfg.pooledZPairs = [ ...
%     "Glu",     "GABA"; ...
%     "Glu+Gln", "GABA"];

% Ranked table and exports
plotCfg.doRankedTable = true;
plotCfg.exportRankedTable = true;
plotCfg.exportCRLBQualityTable = true;

% Saving figures
% plotCfg.saveFigures = true;
% plotCfg.outputDir = fullfile(pwd, "TemporalCorrelationPlots");

% Run plotting pipeline
plotOutputs = PlotTemporalCorrelationDiagnostics( ...
    covOutputs, ...
    deGraafOutputs, ...
    plotCfg);


%% Plot CRLB histograms per metabolite

crlbCfg = struct();

% Use the same sum-preferred panel logic.
crlbCfg.metabolites = "all";

% These are the sums we prefer to keep.
crlbCfg.sumMetabolites = ["GPC+PCh", "NAA+NAAG", "Cr+PCr", "Glu+Gln"];

% If true: if NAA+NAAG exists, remove NAA and NAAG, etc.
crlbCfg.useSumPreferredMetabolites = true;

% Histogram settings.
crlbCfg.binWidth = 5;
crlbCfg.crlbLimitLine = 100;

% Output.
crlbCfg.saveFigures = true;
crlbCfg.outputDir = fullfile(pwd, "CRLB_Histograms");

crlbOutputs = PlotCRLBHistogramsFromDeGraafOutputs(deGraafOutputs, crlbCfg);

disp("CRLB summary table:")
disp(crlbOutputs.summaryTable)


%% ================= ANALYSIS FILTERING =================

filterCfg = struct();
filterCfg.patientIDs = "all";
filterCfg.division = 1;
filterCfg.metabolites = "all";

% Preserve the current sum-preferred policy: sums remain and components are
% removed before global CRLB-majority screening.
filterCfg.useSumPreferredFilter = true;
filterCfg.sumMetabolites = ["GPC+PCh", "NAA+NAAG", "Cr+PCr", "Glu+Gln"];

% Pool finite Division-1 CRLB values across common patients and parts.
filterCfg.useCRLBMajorityFilter = true;
filterCfg.crlbMajorityThreshold = 100;

% Canonical empirical observation handling.
filterCfg.ignoreZeros = true;
filterCfg.pairwiseMinValidParts = 10;
filterCfg.pairwiseMinPatients = 3;

% Preserve the temporal test's existing, distinct eligibility threshold.
filterCfg.temporalMinValidParts = 8;
filterCfg.temporalMinPatients = 3;
filterCfg.temporalUseGlobalMetabolites = false;
filterCfg.temporalCRLBThreshold = 100;
filterCfg.temporalRequiredGoodFraction = 0.01;
filterCfg.prepareTemporalCircularShift = true;

% Preserve each Wishart example's requested panel before applying the same
% sum-preferred and global CRLB-majority policy.
filterCfg.wishartMinValidParts = 30;
filterCfg.wishartViews.modeA = struct('metabolites', "all", 'minValidParts', 30);
filterCfg.wishartViews.modeB = struct( ...
    'metabolites', ["NAA", "Cr", "PCr", "Glu"], 'minValidParts', 30);
filterCfg.wishartViews.modeC = struct('metabolites', "all", 'minValidParts', 30);

[analysisData, filterReport] = ApplyAnalysisFilters( ...
    covOutputs, deGraafOutputs, filterCfg);


%% ================= STATISTICAL TEST CONFIGURATION =================

statsCfg = struct();
statsCfg.pairwise = struct();
statsCfg.temporal = struct();
statsCfg.wishart = struct();

% Pairwise one-sample t-test significance threshold. Benjamini-Hochberg
% FDR correction remains enabled by the validated implementation.
statsCfg.pairwise.alpha = 0.05;
statsCfg.pairwise.useFDR = true;

% Preserve the existing pairwise exports and output location.
statsCfg.pairwise.exportResults = true;
statsCfg.pairwise.outputDir = fullfile( ...
    pwd, "PairwiseEmpiricalVsModelResults");

% Wishart covariance LRT mathematical and reporting parameters. Patient and
% metabolite-panel eligibility remains in filterCfg.wishartViews.
statsCfg.wishart.minMetabolites = 2;
statsCfg.wishart.alpha = 0.05;
statsCfg.wishart.applyNumericalRidge = true;
statsCfg.wishart.ridgeScale = 1e-8;
statsCfg.wishart.maxRidgeSteps = 10;
statsCfg.wishart.runOnlyCommonPanelPatients = true;
statsCfg.wishart.exportResults = true;
statsCfg.wishart.outputDirs = struct();
statsCfg.wishart.outputDirs.modeA = fullfile( ...
    pwd, "WishartLRTResults_ModeA_PerPatientLargestValid");
statsCfg.wishart.outputDirs.modeB = fullfile( ...
    pwd, "WishartLRTResults_ModeB_PredefinedPanel");
statsCfg.wishart.outputDirs.modeC = fullfile( ...
    pwd, "WishartLRTResults_ModeC_LargestCommon");


%% Configure temporal statistical tests

statsCfg.temporal.doFisherGroupTest = true;
statsCfg.temporal.doCircularShiftPermutation = true;
statsCfg.temporal.nGroupPermutations = 5000;
statsCfg.temporal.rngSeed = 1;
statsCfg.temporal.useFDR = true;
statsCfg.temporal.exportTables = true;
statsCfg.temporal.outputDir = fullfile(pwd, "TemporalCorrelationStats");

metabs = analysisData.temporal.metabolites;

if isempty(metabs) || all(metabs == "")
    metabs = string(covOutputs.group.meanCorrTable.Properties.VariableNames);
end

nMetabs = numel(metabs);

pairA = strings(0, 1);
pairB = strings(0, 1);

for iMetab = 1:nMetabs
    for jMetab = iMetab+1:nMetabs
        pairA(end+1, 1) = metabs(iMetab);
        pairB(end+1, 1) = metabs(jMetab);
    end
end

statsCfg.temporal.permutationPairs = [pairA, pairB];

temporal = ProjectStatistics.TemporalCorrelation( ...
    analysisData.temporal, statsCfg.temporal);

if ~isempty(temporal.permutationResults)

    S = temporal.permutationResults;

    if ismember("groupCircularShiftPValue", string(S.Properties.VariableNames))

        p = S.groupCircularShiftPValue;
        q = nan(size(p));

        valid = ~isnan(p);
        pValid = p(valid);

        if ~isempty(pValid)

            [pSorted, sortOrder] = sort(pValid, "ascend");
            m = numel(pSorted);

            qSorted = pSorted .* m ./ (1:m)';

            for k = m-1:-1:1
                qSorted(k) = min(qSorted(k), qSorted(k+1));
            end

            qSorted = min(qSorted, 1);

            qValid = nan(size(pValid));
            qValid(sortOrder) = qSorted;

            q(valid) = qValid;
        end

        S.qValueCircularShift_FDR = q;

        if ismember("qValueCircularShift_FDR", string(S.Properties.VariableNames))
            S = sortrows(S, "qValueCircularShift_FDR", "ascend");
        end

        temporal.permutationResults = S;
    end
end

outputDir = fullfile(pwd, "TemporalCorrelationStats_Simplified_AllPairs");

if ~isfolder(outputDir)
    mkdir(outputDir);
end

T = temporal.summaryTable;

wantedColsMain = [ ...
    "metaboliteA", ...
    "metaboliteB", ...
    "groupMeanR", ...
    "pValue", ...
    "qValue_FDR", ...
    "nPatientsUsed", ...
    "nPositivePatients", ...
    "nNegativePatients", ...
    "signConsistency", ...
    "CRLB_pair_status", ...
    "absLCModelCorr"];

availableColsMain = string(T.Properties.VariableNames);
colsToKeepMain = wantedColsMain(ismember(wantedColsMain, availableColsMain));

T_simple = T(:, colsToKeepMain);

if ismember("qValue_FDR", string(T_simple.Properties.VariableNames))
    T_simple = sortrows(T_simple, "qValue_FDR", "ascend");
end

writetable(T_simple, ...
    fullfile(outputDir, "Main_FisherZ_Results_Simplified.csv"));

if ismember("qValue_FDR", availableColsMain) && ...
        ismember("CRLB_pair_status", availableColsMain) && ...
        ismember("signConsistency", availableColsMain)

    T_strong = T( ...
        T.qValue_FDR < 0.05 & ...
        T.CRLB_pair_status == "PASS" & ...
        T.signConsistency >= 0.75, :);

    T_strong = T_strong(:, colsToKeepMain);

    if height(T_strong) > 0 && ...
            ismember("qValue_FDR", string(T_strong.Properties.VariableNames))
        T_strong = sortrows(T_strong, "qValue_FDR", "ascend");
    end

    writetable(T_strong, ...
        fullfile(outputDir, "Strong_Candidates_q_FDR_CRLB_PASS.csv"));
end

if ismember("pValue", availableColsMain) && ...
        ismember("signConsistency", availableColsMain)

    T_exploratory = T( ...
        T.pValue < 0.05 & ...
        T.signConsistency >= 0.75, :);

    T_exploratory = T_exploratory(:, colsToKeepMain);

    if height(T_exploratory) > 0 && ...
            ismember("pValue", string(T_exploratory.Properties.VariableNames))
        T_exploratory = sortrows(T_exploratory, "pValue", "ascend");
    end

    writetable(T_exploratory, ...
        fullfile(outputDir, "Exploratory_Candidates_p_uncorrected.csv"));
end

if ~isempty(temporal.permutationResults)

    S = temporal.permutationResults;

    wantedColsShift = [ ...
        "metaboliteA", ...
        "metaboliteB", ...
        "observedGroupR", ...
        "groupCircularShiftPValue", ...
        "qValueCircularShift_FDR", ...
        "nPatientsUsed"];

    availableColsShift = string(S.Properties.VariableNames);
    colsToKeepShift = wantedColsShift(ismember(wantedColsShift, availableColsShift));

    S_simple = S(:, colsToKeepShift);

    if ismember("qValueCircularShift_FDR", string(S_simple.Properties.VariableNames))
        S_simple = sortrows(S_simple, "qValueCircularShift_FDR", "ascend");
    elseif ismember("groupCircularShiftPValue", string(S_simple.Properties.VariableNames))
        S_simple = sortrows(S_simple, "groupCircularShiftPValue", "ascend");
    end

    writetable(S_simple, ...
        fullfile(outputDir, "Circular_Shift_All_Pairs_Simplified.csv"));
end

if ~isempty(temporal.patientTable)

    P = temporal.patientTable;

    wantedColsPatient = [ ...
        "metaboliteA", ...
        "metaboliteB", ...
        "patientID", ...
        "rValue", ...
        "nValidParts"];

    availableColsPatient = string(P.Properties.VariableNames);
    colsToKeepPatient = wantedColsPatient(ismember(wantedColsPatient, availableColsPatient));

    P_simple = P(:, colsToKeepPatient);

    writetable(P_simple, ...
        fullfile(outputDir, "Patient_Level_Correlations_Simplified.csv"));
end

if ~isempty(temporal.crlbQualityTable)

    Q = temporal.crlbQualityTable;

    wantedColsCRLB = [ ...
        "metabolite", ...
        "nCRLBUnder100", ...
        "nInstances", ...
        "fractionCRLBUnder100", ...
        "fails90PercentRule"];

    availableColsCRLB = string(Q.Properties.VariableNames);
    colsToKeepCRLB = wantedColsCRLB(ismember(wantedColsCRLB, availableColsCRLB));

    Q_simple = Q(:, colsToKeepCRLB);

    if ismember("fractionCRLBUnder100", string(Q_simple.Properties.VariableNames))
        Q_simple = sortrows(Q_simple, "fractionCRLBUnder100", "ascend");
    end

    writetable(Q_simple, ...
        fullfile(outputDir, "CRLB_Reliability_Simplified.csv"));
end

fprintf("\nFinished temporal correlation tests.\n");
fprintf("Number of metabolite pairs tested with Fisher-z: %d\n", height(temporal.summaryTable));

if ~isempty(temporal.permutationResults)
    fprintf("Number of metabolite pairs tested with circular shift: %d\n", height(temporal.permutationResults));
end

fprintf("Saved simplified CSV files to:\n%s\n", outputDir);



%% Wishart covariance LRT

wishart = ProjectStatistics.WishartCovarianceLRT( ...
    analysisData.wishart, statsCfg.wishart);

%% Mode A: per-patient largest valid subset
% Each patient gets the largest valid non-sum metabolite subset available for that patient.
disp("Mode A: per-patient largest valid subset")
disp(wishart.modes.modeA.patientTable)
disp(wishart.modes.modeA.summaryTable)

%% Mode B: predefined fixed metabolite panel
% The exact same predefined metabolite panel is requested for every patient.
disp("Mode B: predefined fixed metabolite panel")
disp("Fixed panel used:")
disp(wishart.modes.modeB.diagnostics.preparedMetabolites')
disp(wishart.modes.modeB.patientTable)
disp(wishart.modes.modeB.summaryTable)

%% Mode C: largest common valid subset
% First find the largest valid panel for each patient, then take the common
% intersection, then run the LRT using that same common panel for every eligible patient.
disp("Mode C: largest common valid subset")
disp("Common panel used:")
disp(wishart.modes.modeC.commonPanel')
disp(wishart.modes.modeC.panelDiscoveryTable)
disp(wishart.modes.modeC.patientTable)
disp(wishart.modes.modeC.summaryTable)


%% Pairwise empirical-vs-LCModel/de Graaf correlation test
% Option A: one-sample t-test across patients on Fisher-z differences

% Run test.
pairwise = ProjectStatistics.PairwiseEmpiricalVsModel( ...
    analysisData.pairwise, statsCfg.pairwise);

%% Display results

disp("Pairwise empirical-vs-LCModel/deGraaf summary:")
disp(pairwise.summaryTable)

disp("Patient-level pairwise results:")
disp(pairwise.patientTable)

%% Show only FDR-significant pairs

sigPairs = pairwise.summaryTable(pairwise.summaryTable.rejectH0_FDR == true, :);

disp("FDR-significant empirical-vs-model pairs:")
disp(sigPairs)
