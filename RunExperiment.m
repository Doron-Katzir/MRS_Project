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
% cfg.input.singlePatientID = "P11";
% cfg.input.singleFileName = "meas_MID00090_FID32072_eja_svs_slaser_TE_80_r0.dat";
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
outputs = Load_subsets_multi_patient_input(cfg);


%% Check metab covariance

cfg.covariance.loadMode = "allSubfolders";
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
cfg.degraaf.maskInvalidCRLB = true;
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
plotCfg.requiredGoodFraction = 0.90;

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


%% Configure statistical tests

testCfg = struct();

testCfg.minValidParts = 8;
testCfg.minPatientsForGroupTest = 3;
testCfg.nGroupPermutations = 5000;
testCfg.randomSeed = 1;
testCfg.crlbThreshold = 100;
testCfg.requiredGoodFraction = 0.90;

metabs = string(covOutputs.group.meanCorrTable.Properties.RowNames);

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

testCfg.permutationPairs = [pairA, pairB];

testOutputs = TestTemporalMetaboliteCorrelations( ...
    covOutputs, ...
    deGraafOutputs, ...
    testCfg);

if isfield(testOutputs, "circularShiftTable") && ...
        ~isempty(testOutputs.circularShiftTable)

    S = testOutputs.circularShiftTable;

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

        testOutputs.circularShiftTable = S;
    end
end

outputDir = fullfile(pwd, "TemporalCorrelationStats_Simplified_AllPairs");

if ~isfolder(outputDir)
    mkdir(outputDir);
end

T = testOutputs.groupFisherTable;

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

if isfield(testOutputs, "circularShiftTable") && ...
        ~isempty(testOutputs.circularShiftTable)

    S = testOutputs.circularShiftTable;

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

if isfield(testOutputs, "patientCorrelationTable") && ...
        ~isempty(testOutputs.patientCorrelationTable)

    P = testOutputs.patientCorrelationTable;

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

if isfield(testOutputs, "groupCRLBQualityTable") && ...
        ~isempty(testOutputs.groupCRLBQualityTable)

    Q = testOutputs.groupCRLBQualityTable;

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
fprintf("Number of metabolite pairs tested with Fisher-z: %d\n", height(testOutputs.groupFisherTable));

if isfield(testOutputs, "circularShiftTable") && ...
        ~isempty(testOutputs.circularShiftTable)
    fprintf("Number of metabolite pairs tested with circular shift: %d\n", height(testOutputs.circularShiftTable));
end

fprintf("Saved simplified CSV files to:\n%s\n", outputDir);



%% Wishart / covariance LRT examples
% Put one of these blocks in RunExperiment.m after:
%   covOutputs = MetabCovarianceByPatient(cfg);
%   deGraafOutputs = DeGraafAmplitudeCorrelationByPatient(cfg);

% Shared settings
baseLrtCfg = struct();
baseLrtCfg.patientIDs = "all";
baseLrtCfg.excludeSumMetabolites = true;
baseLrtCfg.sumMetabolites = ["GPC+PCh", "NAA+NAAG", "Cr+PCr", "Glu+Gln"];
baseLrtCfg.minValidParts = 10;
baseLrtCfg.minMetabolites = 2;
baseLrtCfg.alpha = 0.05;
baseLrtCfg.applyNumericalRidge = true;
baseLrtCfg.ridgeScale = 1e-8;
baseLrtCfg.exportResults = true;

%% Mode A: per-patient largest valid subset
% Each patient gets the largest valid non-sum metabolite subset available for that patient.
lrtCfg = baseLrtCfg;
lrtCfg.metaboliteSelectionMode = "perPatientLargestValid";  % aliases: "perPatient" or "a"
lrtCfg.metabolites = "all";
lrtCfg.outputDir = fullfile(pwd, "WishartLRTResults_ModeA_PerPatientLargestValid");

lrtOutputs_ModeA = TestWishartCovarianceLRT(covOutputs, deGraafOutputs, lrtCfg);

disp("Mode A: per-patient largest valid subset")
disp(lrtOutputs_ModeA.patientSummaryTable)
disp(lrtOutputs_ModeA.groupTable)

%% Mode B: predefined fixed metabolite panel
% The exact same predefined metabolite panel is requested for every patient.
lrtCfg = baseLrtCfg;
lrtCfg.metaboliteSelectionMode = "fixed";  % aliases: "predefined" or "b"
lrtCfg.metabolites = ["NAA", "Cr", "PCr", "Glu"];
lrtCfg.outputDir = fullfile(pwd, "WishartLRTResults_ModeB_PredefinedPanel");

lrtOutputs_ModeB = TestWishartCovarianceLRT(covOutputs, deGraafOutputs, lrtCfg);

disp("Mode B: predefined fixed metabolite panel")
disp("Fixed panel used:")
disp(lrtCfg.metabolites')
disp(lrtOutputs_ModeB.patientSummaryTable)
disp(lrtOutputs_ModeB.groupTable)

%% Mode C: largest common valid subset
% First find the largest valid panel for each patient, then take the common
% intersection, then run the LRT using that same common panel for every eligible patient.
lrtCfg = baseLrtCfg;
lrtCfg.metaboliteSelectionMode = "largestCommon";  % aliases: "common" or "c"
lrtCfg.metabolites = "all";
lrtCfg.runOnlyCommonPanelPatients = true;
lrtCfg.outputDir = fullfile(pwd, "WishartLRTResults_ModeC_LargestCommon");

lrtOutputs_ModeC = TestWishartCovarianceLRT(covOutputs, deGraafOutputs, lrtCfg);

disp("Mode C: largest common valid subset")
disp("Common panel used:")
disp(lrtOutputs_ModeC.commonPanel')
disp(lrtOutputs_ModeC.panelDiscoveryTable)
disp(lrtOutputs_ModeC.patientSummaryTable)
disp(lrtOutputs_ModeC.groupTable)


%% Pairwise empirical-vs-LCModel/de Graaf correlation test
% Option A: one-sample t-test across patients on Fisher-z differences

pairCfg = struct();

% Use all patients present in both outputs.
pairCfg.patientIDs = "all";

% Test all available non-sum metabolites.
pairCfg.metabolites = "all";

% Or manually use a cleaner fixed panel:
% pairCfg.metabolites = ["NAA", "Cr", "PCr", "Glu"];

% Exclude summed metabolites.
pairCfg.excludeSumMetabolites = true;
pairCfg.sumMetabolites = ["GPC+PCh", "NAA+NAAG", "Cr+PCr", "Glu+Gln"];

% Pairwise valid repeated parts.
pairCfg.ignoreZeros = true;
pairCfg.minValidParts = 10;

% Minimum patients required for a group-level test.
pairCfg.minPatientsForGroupTest = 3;

% Significance threshold.
pairCfg.alpha = 0.05;

% Export results.
pairCfg.exportResults = true;
pairCfg.outputDir = fullfile(pwd, "PairwiseEmpiricalVsModelResults");

% Run test.
pairOutputs = TestPairwiseEmpiricalVsModelCorrelation(covOutputs, deGraafOutputs, pairCfg);

%% Display results

disp("Pairwise empirical-vs-LCModel/deGraaf summary:")
disp(pairOutputs.pairSummaryTable)

disp("Patient-level pairwise results:")
disp(pairOutputs.patientPairTable)

%% Show only FDR-significant pairs

sigPairs = pairOutputs.pairSummaryTable(pairOutputs.pairSummaryTable.rejectH0_FDR == true, :);

disp("FDR-significant empirical-vs-model pairs:")
disp(sigPairs)
