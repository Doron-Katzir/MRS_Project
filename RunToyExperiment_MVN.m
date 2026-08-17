%% RunToyExperiment_MVN_Clean.m
%
% Clean toy validation script for the MRS temporal-correlation pipeline.
%
% Purpose:
%   Generate fake repeated metabolite estimates that follow the assumptions
%   behind the Wishart covariance LRT and Fisher-z pairwise correlation tests,
%   then run the SAME CORE STATISTICAL FUNCTIONS used in the real analysis:
%
%       1) TestWishartCovarianceLRT
%       2) TestPairwiseEmpiricalVsModelCorrelation
%
% The toy script is intentionally much cleaner than RunExperiment.m.
% It does NOT run all exploratory diagnostics. It only creates the outputs
% needed to answer:
%
%   "Does the statistical pipeline behave as expected on data where we know
%    the truth?"
%
% Main outputs:
%   Toy_MVN_Clean_Output/
%       Toy_Run_Summary.csv
%       Toy_Matrix_LRT_Summary.csv
%       Toy_Pairwise_Summary.csv
%       Toy_Truth_Table.csv
%       Figure_1_Toy_Truth_Matrices.png
%       Figure_2_Null_Observed_vs_Model.png
%       Figure_3_Alternative_Observed_vs_Model.png
%
% Optional debug files are disabled by default.
%
% Assumptions:
%   - TestWishartCovarianceLRT.m is on the MATLAB path.
%   - TestPairwiseEmpiricalVsModelCorrelation.m is on the MATLAB path.

clear;
clc;
close all;

%% ------------------------------------------------------------------------
% User settings
% -------------------------------------------------------------------------

toyCfg = struct();

toyCfg.randomSeed = 1;

toyCfg.nPatients = 53;
toyCfg.nParts = 36;

% Keep the toy panel small and clean.
% This matches the common-panel size in the real analysis.
toyCfg.metabolites = ["NAA", "Cr", "PCr", "Glu"];

toyCfg.alpha = 0.05;

toyCfg.outputDir = fullfile(pwd, "Toy_MVN_Clean_Output");

% Keep the main output clean.
% Change to true only if you want patient-level debug tables/logs.
toyCfg.saveDebugFiles = false;

% Delete previous toy-output folder so that results do not accumulate.
toyCfg.deleteOldOutputDir = true;

%% ------------------------------------------------------------------------
% Check required functions
% -------------------------------------------------------------------------

requiredFunctions = [
    "TestWishartCovarianceLRT"
    "TestPairwiseEmpiricalVsModelCorrelation"
    ];

for k = 1:numel(requiredFunctions)
    if exist(requiredFunctions(k), "file") ~= 2
        error("Missing required function on MATLAB path: %s.m", requiredFunctions(k));
    end
end

%% ------------------------------------------------------------------------
% Prepare output folder
% -------------------------------------------------------------------------

if toyCfg.deleteOldOutputDir && isfolder(toyCfg.outputDir)
    rmdir(toyCfg.outputDir, "s");
end

if ~isfolder(toyCfg.outputDir)
    mkdir(toyCfg.outputDir);
end

debugDir = fullfile(toyCfg.outputDir, "Debug_Optional");
if toyCfg.saveDebugFiles && ~isfolder(debugDir)
    mkdir(debugDir);
end

%% ------------------------------------------------------------------------
% Define the fake LCModel/de Graaf model covariance
% -------------------------------------------------------------------------

metabs = toyCfg.metabolites(:);
p = numel(metabs);

% Fake model standard deviations.
% These represent model uncertainty scale for the fake LCModel/de Graaf model.
modelSD = [2.0, 1.0, 0.9, 1.5];

% Fake LCModel/de Graaf model correlation matrix.
% This is the model that the tests compare against.
R_model = [
     1.00,   0.10,  -0.10,   0.05
     0.10,   1.00,  -0.65,   0.10
    -0.10,  -0.65,   1.00,  -0.05
     0.05,   0.10,  -0.05,   1.00
    ];

Sigma_model = diag(modelSD) * R_model * diag(modelSD);

AssertPositiveDefinite(R_model, "R_model");
AssertPositiveDefinite(Sigma_model, "Sigma_model");

%% ------------------------------------------------------------------------
% Define two toy scenarios
% -------------------------------------------------------------------------

% Scenario 1: null
% The true covariance equals the model covariance.
R_null = R_model;
Sigma_null = Sigma_model;

% Scenario 2: alternative
% The true covariance intentionally differs from the model covariance.
% We plant two known differences:
%   NAA-Glu: model 0.05 -> true 0.65
%   Cr-PCr : model -0.65 -> true -0.90
R_alt = R_model;

idxNAA = FindMetabIndex(metabs, "NAA");
idxCr  = FindMetabIndex(metabs, "Cr");
idxPCr = FindMetabIndex(metabs, "PCr");
idxGlu = FindMetabIndex(metabs, "Glu");

R_alt(idxNAA, idxGlu) = 0.65;
R_alt(idxGlu, idxNAA) = 0.65;

R_alt(idxCr, idxPCr) = -0.90;
R_alt(idxPCr, idxCr) = -0.90;

Sigma_alt = diag(modelSD) * R_alt * diag(modelSD);

AssertPositiveDefinite(R_alt, "R_alt");
AssertPositiveDefinite(Sigma_alt, "Sigma_alt");

scenarios = struct([]);

scenarios(1).name = "null";
scenarios(1).description = "Sigma_true equals Sigma_model";
scenarios(1).SigmaTrue = Sigma_null;
scenarios(1).RTrue = R_null;
scenarios(1).expectedMatrixReject = false;
scenarios(1).expectedChangedPairs = false(p, p);

scenarios(2).name = "alternative";
scenarios(2).description = "Sigma_true differs from Sigma_model in selected pairs";
scenarios(2).SigmaTrue = Sigma_alt;
scenarios(2).RTrue = R_alt;
scenarios(2).expectedMatrixReject = true;
scenarios(2).expectedChangedPairs = abs(R_alt - R_model) > 1e-12;

%% ------------------------------------------------------------------------
% Plot the planted truth
% -------------------------------------------------------------------------

PlotToyTruthMatrices( ...
    metabs, ...
    R_model, ...
    R_null, ...
    R_alt, ...
    fullfile(toyCfg.outputDir, "Figure_1_Toy_Truth_Matrices.png"));

%% ------------------------------------------------------------------------
% Run both scenarios
% -------------------------------------------------------------------------

matrixSummaryAll = table();
pairwiseSummaryAll = table();
runSummaryAll = table();
truthTableAll = table();

for sIdx = 1:numel(scenarios)

    scenario = scenarios(sIdx);

    fprintf("\n============================================================\n");
    fprintf("Running toy scenario: %s\n", scenario.name);
    fprintf("%s\n", scenario.description);
    fprintf("============================================================\n");

    scenarioOutputDir = fullfile(toyCfg.outputDir, "Scenario_" + scenario.name);
    if ~isfolder(scenarioOutputDir)
        mkdir(scenarioOutputDir);
    end

    rng(toyCfg.randomSeed + sIdx - 1);

    [covOutputs, deGraafOutputs] = BuildToyPipelineStructures( ...
        toyCfg.nPatients, ...
        toyCfg.nParts, ...
        metabs, ...
        scenario.SigmaTrue, ...
        Sigma_model);

    truthTable = BuildTruthTable( ...
        scenario.name, ...
        metabs, ...
        R_model, ...
        scenario.RTrue, ...
        scenario.expectedChangedPairs);

    [matrixSummary, lrtOutputs, lrtLogText] = RunCleanWishartLRT( ...
        covOutputs, ...
        deGraafOutputs, ...
        toyCfg, ...
        scenario);

    [pairwiseSummary, pairOutputs, pairLogText] = RunCleanPairwiseTest( ...
        covOutputs, ...
        deGraafOutputs, ...
        toyCfg, ...
        scenario, ...
        truthTable);

    observedR = table2array(covOutputs.group.meanCorrTable);

    if scenario.name == "null"
        figFile = fullfile(toyCfg.outputDir, "Figure_2_Null_Observed_vs_Model.png");
    else
        figFile = fullfile(toyCfg.outputDir, "Figure_3_Alternative_Observed_vs_Model.png");
    end

    PlotObservedVsModelMatrices( ...
        metabs, ...
        R_model, ...
        observedR, ...
        scenario.RTrue, ...
        scenario.name, ...
        figFile);

    runSummary = BuildScenarioRunSummary( ...
        scenario, ...
        matrixSummary, ...
        pairwiseSummary, ...
        toyCfg);

    matrixSummaryAll = [matrixSummaryAll; matrixSummary]; %#ok<AGROW>
    pairwiseSummaryAll = [pairwiseSummaryAll; pairwiseSummary]; %#ok<AGROW>
    runSummaryAll = [runSummaryAll; runSummary]; %#ok<AGROW>
    truthTableAll = [truthTableAll; truthTable]; %#ok<AGROW>

    if toyCfg.saveDebugFiles
        writetable(lrtOutputs.patientSummaryTable, ...
            fullfile(debugDir, "Debug_" + scenario.name + "_Patient_Level_LRT.csv"));

        writetable(pairOutputs.patientPairTable, ...
            fullfile(debugDir, "Debug_" + scenario.name + "_Patient_Level_Pairwise.csv"));

        WriteTextFile(fullfile(debugDir, "Debug_" + scenario.name + "_LRT_log.txt"), lrtLogText);
        WriteTextFile(fullfile(debugDir, "Debug_" + scenario.name + "_Pairwise_log.txt"), pairLogText);
    end
end

%% ------------------------------------------------------------------------
% Export only clean summary files
% -------------------------------------------------------------------------

writetable(runSummaryAll, ...
    fullfile(toyCfg.outputDir, "Toy_Run_Summary.csv"));

writetable(matrixSummaryAll, ...
    fullfile(toyCfg.outputDir, "Toy_Matrix_LRT_Summary.csv"));

writetable(pairwiseSummaryAll, ...
    fullfile(toyCfg.outputDir, "Toy_Pairwise_Summary.csv"));

writetable(truthTableAll, ...
    fullfile(toyCfg.outputDir, "Toy_Truth_Table.csv"));

%% ------------------------------------------------------------------------
% Print compact console summary
% -------------------------------------------------------------------------

fprintf("\n============================================================\n");
fprintf("TOY VALIDATION COMPLETE\n");
fprintf("============================================================\n");

fprintf("\nRun summary:\n");
disp(runSummaryAll)

fprintf("\nMatrix LRT summary:\n");
disp(matrixSummaryAll)

fprintf("\nPairwise summary, sorted by scenario and q-value:\n");

pairDisplay = pairwiseSummaryAll;
if ismember("qValue_FDR", string(pairDisplay.Properties.VariableNames))
    pairDisplay = sortrows(pairDisplay, ["scenario", "qValue_FDR"], ["ascend", "ascend"]);
end

disp(pairDisplay)

fprintf("\nClean toy outputs saved to:\n%s\n", toyCfg.outputDir);

fprintf("\nGenerated files:\n");
fprintf("  Toy_Run_Summary.csv\n");
fprintf("  Toy_Matrix_LRT_Summary.csv\n");
fprintf("  Toy_Pairwise_Summary.csv\n");
fprintf("  Toy_Truth_Table.csv\n");
fprintf("  Figure_1_Toy_Truth_Matrices.png\n");
fprintf("  Figure_2_Null_Observed_vs_Model.png\n");
fprintf("  Figure_3_Alternative_Observed_vs_Model.png\n");

if toyCfg.saveDebugFiles
    fprintf("  Debug_Optional/...\n");
end

%% =========================================================================
% Local functions
% =========================================================================

function [covOutputs, deGraafOutputs] = BuildToyPipelineStructures( ...
    nPatients, nParts, metabs, SigmaTrue, SigmaModel)

    metabs = string(metabs(:));
    nMetabs = numel(metabs);

    RModel = CovToCorr(SigmaModel);

    covOutputs = struct();
    deGraafOutputs = struct();

    covOutputs.patientResultsByID = struct();
    deGraafOutputs.patientResultsByID = struct();

    patientIDs = strings(nPatients, 1);

    empCovStack = nan(nMetabs, nMetabs, nPatients);
    empCorrStack = nan(nMetabs, nMetabs, nPatients);

    modelCovStack = nan(nMetabs, nMetabs, nPatients);
    modelCorrStack = nan(nMetabs, nMetabs, nPatients);

    baseMu = BuildBaseMeans(metabs);

    Ltrue = chol(SigmaTrue, "lower");

    modelCovTable = MatrixToSquareTable(SigmaModel, metabs);
    modelCorrTable = MatrixToSquareTable(RModel, metabs);

    for pIdx = 1:nPatients

        patientID = "P" + sprintf("%02d", pIdx);
        patientField = matlab.lang.makeValidName(char(patientID));
        patientIDs(pIdx) = patientID;

        % Small patient-specific baseline shift.
        % This affects means, not the within-patient covariance structure.
        patientMu = baseMu + 0.25 * randn(1, nMetabs);

        Z = randn(nParts, nMetabs);
        X = Z * Ltrue.' + patientMu;

        partTable = array2table(X, "VariableNames", cellstr(metabs));
        partTable.part = (1:nParts).';
        partTable = movevars(partTable, "part", "Before", 1);

        Semp = cov(X, 0);
        Remp = corrcoef(X);

        covOutputs.patientResultsByID.(patientField).patientID = patientID;
        covOutputs.patientResultsByID.(patientField).partTable = partTable;
        covOutputs.patientResultsByID.(patientField).covTable = MatrixToSquareTable(Semp, metabs);
        covOutputs.patientResultsByID.(patientField).corrTable = MatrixToSquareTable(Remp, metabs);
        covOutputs.patientResultsByID.(patientField).nPairTable = MatrixToSquareTable(nParts * ones(nMetabs), metabs);

        deGraafOutputs.patientResultsByID.(patientField).patientID = patientID;
        deGraafOutputs.patientResultsByID.(patientField).meanAmplitudeCovTable = modelCovTable;
        deGraafOutputs.patientResultsByID.(patientField).meanAmplitudeCorrTable = modelCorrTable;
        deGraafOutputs.patientResultsByID.(patientField).crlbSummaryTable = BuildToyCRLBSummary(metabs, nParts);

        empCovStack(:, :, pIdx) = Semp;
        empCorrStack(:, :, pIdx) = Remp;
        modelCovStack(:, :, pIdx) = SigmaModel;
        modelCorrStack(:, :, pIdx) = RModel;
    end

    covOutputs.patientIDs = patientIDs;
    covOutputs.group.meanCovTable = MatrixToSquareTable(mean(empCovStack, 3), metabs);
    covOutputs.group.meanCorrTable = MatrixToSquareTable(FisherMeanCorrStack(empCorrStack), metabs);
    covOutputs.group.meanAbsCorrTable = MatrixToSquareTable(mean(abs(empCorrStack), 3), metabs);
    covOutputs.group.nPatientsCovTable = MatrixToSquareTable(nPatients * ones(nMetabs), metabs);
    covOutputs.group.nPatientsCorrTable = MatrixToSquareTable(nPatients * ones(nMetabs), metabs);

    deGraafOutputs.patientIDs = patientIDs;
    deGraafOutputs.group.meanAmplitudeCovTable = MatrixToSquareTable(mean(modelCovStack, 3), metabs);
    deGraafOutputs.group.meanAmplitudeCorrTable = MatrixToSquareTable(FisherMeanCorrStack(modelCorrStack), metabs);
    deGraafOutputs.group.meanAbsAmplitudeCorrTable = MatrixToSquareTable(mean(abs(modelCorrStack), 3), metabs);
    deGraafOutputs.group.nPatientsCovTable = MatrixToSquareTable(nPatients * ones(nMetabs), metabs);
    deGraafOutputs.group.nPatientsCorrTable = MatrixToSquareTable(nPatients * ones(nMetabs), metabs);
    deGraafOutputs.group.crlbQualityTable = BuildToyCRLBQuality(metabs, nPatients, nParts);
end

function [matrixSummary, lrtOutputs, capturedText] = RunCleanWishartLRT( ...
    covOutputs, deGraafOutputs, toyCfg, scenario)

    lrtCfg = struct();

    lrtCfg.patientIDs = "all";

    % For toy validation, use one fixed panel.
    % This keeps the result interpretable.
    lrtCfg.metaboliteSelectionMode = "fixed";
    lrtCfg.metabolites = toyCfg.metabolites;

    lrtCfg.excludeSumMetabolites = true;
    lrtCfg.sumMetabolites = ["GPC+PCh", "NAA+NAAG", "Cr+PCr", "Glu+Gln"];

    lrtCfg.minValidParts = 10;
    lrtCfg.minMetabolites = 2;

    lrtCfg.alpha = toyCfg.alpha;

    lrtCfg.applyNumericalRidge = true;
    lrtCfg.ridgeScale = 1e-8;

    % We export our own clean summary files from this toy script.
    lrtCfg.exportResults = false;
    lrtCfg.outputDir = fullfile(toyCfg.outputDir, "internal_lrt_not_used");

    [capturedText, lrtOutputs] = CaptureFunctionCall( ...
        @() TestWishartCovarianceLRT(covOutputs, deGraafOutputs, lrtCfg));

    G = lrtOutputs.groupTable;

    nPatientsUsed = ExtractNumeric(G, "nPatientsUsed");
    Tgroup = ExtractNumeric(G, "T_group");
    dfGroup = ExtractNumeric(G, "df_group");
    pValueGroup = ExtractNumeric(G, "pValue_group");
    rejectGroup = ExtractLogical(G, "rejectH0_group");

    expectedReject = logical(scenario.expectedMatrixReject);
    behavior = BehaviorLabel(rejectGroup, expectedReject);

    matrixSummary = table( ...
        string(scenario.name), ...
        string(scenario.description), ...
        nPatientsUsed, ...
        toyCfg.nParts, ...
        string(strjoin(cellstr(toyCfg.metabolites(:)), ", ")), ...
        Tgroup, ...
        dfGroup, ...
        pValueGroup, ...
        rejectGroup, ...
        expectedReject, ...
        string(behavior), ...
        'VariableNames', { ...
            'scenario', ...
            'description', ...
            'nPatientsUsed', ...
            'nParts', ...
            'metabolitesUsed', ...
            'T_group', ...
            'df_group', ...
            'pValue_group', ...
            'rejectH0_group', ...
            'expectedReject', ...
            'behavior'});
end

function [pairwiseSummary, pairOutputs, capturedText] = RunCleanPairwiseTest( ...
    covOutputs, deGraafOutputs, toyCfg, scenario, truthTable)

    pairCfg = struct();

    pairCfg.patientIDs = "all";
    pairCfg.metabolites = toyCfg.metabolites;

    pairCfg.excludeSumMetabolites = true;
    pairCfg.sumMetabolites = ["GPC+PCh", "NAA+NAAG", "Cr+PCr", "Glu+Gln"];

    pairCfg.ignoreZeros = true;
    pairCfg.minValidParts = 10;
    pairCfg.minPatientsForGroupTest = 3;

    pairCfg.alpha = toyCfg.alpha;

    % We export our own clean summary files from this toy script.
    pairCfg.exportResults = false;
    pairCfg.outputDir = fullfile(toyCfg.outputDir, "internal_pairwise_not_used");

    [capturedText, pairOutputs] = CaptureFunctionCall( ...
        @() TestPairwiseEmpiricalVsModelCorrelation(covOutputs, deGraafOutputs, pairCfg));

    pairwiseSummary = BuildCleanPairwiseSummary( ...
        scenario.name, ...
        pairOutputs.pairSummaryTable, ...
        truthTable, ...
        toyCfg.alpha);
end

function runSummary = BuildScenarioRunSummary(scenario, matrixSummary, pairwiseSummary, toyCfg)

    expectedChanged = pairwiseSummary.expectedChanged;
    rejected = pairwiseSummary.rejectH0_FDR;

    nExpectedChanged = sum(expectedChanged);
    nDetectedExpectedChanged = sum(expectedChanged & rejected);
    nFalsePositive = sum(~expectedChanged & rejected);
    nMissedChanged = sum(expectedChanged & ~rejected);

    pairwiseBehaviorPass = nFalsePositive == 0 && nMissedChanged == 0;

    runPass = matrixSummary.behavior == "PASS" && pairwiseBehaviorPass;

    if pairwiseBehaviorPass
        pairBehavior = "PASS";
    else
        pairBehavior = "CHECK";
    end

    if runPass
        overallBehavior = "PASS";
    else
        overallBehavior = "CHECK";
    end

    runSummary = table( ...
        string(scenario.name), ...
        string(scenario.description), ...
        toyCfg.nPatients, ...
        toyCfg.nParts, ...
        string(strjoin(cellstr(toyCfg.metabolites(:)), ", ")), ...
        matrixSummary.expectedReject, ...
        matrixSummary.rejectH0_group, ...
        matrixSummary.pValue_group, ...
        matrixSummary.behavior, ...
        nExpectedChanged, ...
        nDetectedExpectedChanged, ...
        nFalsePositive, ...
        nMissedChanged, ...
        string(pairBehavior), ...
        string(overallBehavior), ...
        'VariableNames', { ...
            'scenario', ...
            'description', ...
            'nPatients', ...
            'nParts', ...
            'metabolites', ...
            'matrixExpectedReject', ...
            'matrixObservedReject', ...
            'matrixPValue', ...
            'matrixBehavior', ...
            'nExpectedChangedPairs', ...
            'nDetectedExpectedChangedPairs', ...
            'nFalsePositivePairs', ...
            'nMissedChangedPairs', ...
            'pairwiseBehavior', ...
            'overallBehavior'});
end

function truthTable = BuildTruthTable(scenarioName, metabs, RModel, RTrue, expectedChangedMatrix)

    metabs = string(metabs(:));
    nMetabs = numel(metabs);

    rows = table();

    for i = 1:nMetabs
        for j = i+1:nMetabs

            rModel = RModel(i, j);
            rTrue = RTrue(i, j);
            delta = rTrue - rModel;
            expectedChanged = expectedChangedMatrix(i, j);

            row = table( ...
                string(scenarioName), ...
                metabs(i), ...
                metabs(j), ...
                rModel, ...
                rTrue, ...
                delta, ...
                expectedChanged, ...
                'VariableNames', { ...
                    'scenario', ...
                    'metaboliteA', ...
                    'metaboliteB', ...
                    'r_model', ...
                    'r_true', ...
                    'trueDeltaR', ...
                    'expectedChanged'});

            rows = [rows; row]; %#ok<AGROW>
        end
    end

    truthTable = rows;
end

function cleanTable = BuildCleanPairwiseSummary(scenarioName, pairSummaryTable, truthTable, alpha)

    nRows = height(truthTable);

    scenario = strings(nRows, 1);
    metaboliteA = strings(nRows, 1);
    metaboliteB = strings(nRows, 1);

    r_model = nan(nRows, 1);
    r_true = nan(nRows, 1);
    trueDeltaR = nan(nRows, 1);
    expectedChanged = false(nRows, 1);

    groupEmpiricalR = nan(nRows, 1);
    groupModelR = nan(nRows, 1);
    meanDeltaZ = nan(nRows, 1);
    pValue = nan(nRows, 1);
    qValueFDR = nan(nRows, 1);
    rejectH0FDR = false(nRows, 1);
    signConsistency = nan(nRows, 1);
    behavior = strings(nRows, 1);

    for rIdx = 1:nRows

        metabA = string(truthTable.metaboliteA(rIdx));
        metabB = string(truthTable.metaboliteB(rIdx));

        scenario(rIdx) = string(scenarioName);
        metaboliteA(rIdx) = metabA;
        metaboliteB(rIdx) = metabB;

        r_model(rIdx) = truthTable.r_model(rIdx);
        r_true(rIdx) = truthTable.r_true(rIdx);
        trueDeltaR(rIdx) = truthTable.trueDeltaR(rIdx);
        expectedChanged(rIdx) = truthTable.expectedChanged(rIdx);

        pairRowIdx = FindPairRow(pairSummaryTable, metabA, metabB);

        if isempty(pairRowIdx)
            rejectH0FDR(rIdx) = false;
            behavior(rIdx) = BehaviorLabel(false, expectedChanged(rIdx));
            continue;
        end

        groupEmpiricalR(rIdx) = ExtractNumeric(pairSummaryTable(pairRowIdx, :), "groupEmpiricalR", NaN);
        groupModelR(rIdx) = ExtractNumeric(pairSummaryTable(pairRowIdx, :), "groupModelR", NaN);
        meanDeltaZ(rIdx) = ExtractNumeric(pairSummaryTable(pairRowIdx, :), "meanDeltaZ_empMinusModel", NaN);
        pValue(rIdx) = ExtractNumeric(pairSummaryTable(pairRowIdx, :), "pValue", NaN);
        qValueFDR(rIdx) = ExtractNumeric(pairSummaryTable(pairRowIdx, :), "qValue_FDR", NaN);
        signConsistency(rIdx) = ExtractNumeric(pairSummaryTable(pairRowIdx, :), "signConsistency", NaN);

        if ismember("rejectH0_FDR", string(pairSummaryTable.Properties.VariableNames))
            rejectH0FDR(rIdx) = logical(pairSummaryTable.rejectH0_FDR(pairRowIdx));
        elseif isfinite(qValueFDR(rIdx))
            rejectH0FDR(rIdx) = qValueFDR(rIdx) < alpha;
        else
            rejectH0FDR(rIdx) = false;
        end

        behavior(rIdx) = BehaviorLabel(rejectH0FDR(rIdx), expectedChanged(rIdx));
    end

    cleanTable = table( ...
        scenario, ...
        metaboliteA, ...
        metaboliteB, ...
        r_model, ...
        r_true, ...
        trueDeltaR, ...
        expectedChanged, ...
        groupEmpiricalR, ...
        groupModelR, ...
        meanDeltaZ, ...
        pValue, ...
        qValueFDR, ...
        rejectH0FDR, ...
        signConsistency, ...
        behavior, ...
        'VariableNames', { ...
            'scenario', ...
            'metaboliteA', ...
            'metaboliteB', ...
            'r_model', ...
            'r_true', ...
            'trueDeltaR', ...
            'expectedChanged', ...
            'groupEmpiricalR', ...
            'groupModelR', ...
            'meanDeltaZ_empMinusModel', ...
            'pValue', ...
            'qValue_FDR', ...
            'rejectH0_FDR', ...
            'signConsistency', ...
            'behavior'});
end

function pairRowIdx = FindPairRow(T, metabA, metabB)

    pairRowIdx = [];

    if isempty(T) || height(T) == 0
        return;
    end

    vars = string(T.Properties.VariableNames);
    if ~ismember("metaboliteA", vars) || ~ismember("metaboliteB", vars)
        return;
    end

    A = string(T.metaboliteA);
    B = string(T.metaboliteB);

    pairRowIdx = find((A == metabA & B == metabB) | (A == metabB & B == metabA), 1);
end

function behavior = BehaviorLabel(observedReject, expectedReject)

    observedReject = logical(observedReject);
    expectedReject = logical(expectedReject);

    if observedReject == expectedReject
        behavior = "PASS";
    else
        behavior = "CHECK";
    end
end

function PlotToyTruthMatrices(metabs, RModel, RNull, RAlt, outputFile)

    fig = figure("Name", "Toy truth matrices", "Color", "w", "Position", [100, 100, 1400, 450]);
    tiledlayout(1, 3, "Padding", "compact", "TileSpacing", "compact");

    ax = nexttile;
    PlotHeatmap(ax, RModel, metabs, "Model correlation R_{model}", [-1, 1], "gray");

    ax = nexttile;
    PlotHeatmap(ax, RNull - RModel, metabs, "Null: R_{true} - R_{model}", [-1, 1], "gray");

    ax = nexttile;
    PlotHeatmap(ax, RAlt - RModel, metabs, "Alternative: R_{true} - R_{model}", [-1, 1], "gray");

    sgtitle("Simulation validation: planted truth", ...
        "Interpreter", "none", ...
        "Color", [0 0 0], ...
        "FontWeight", "bold");

    SaveFigure(fig, outputFile);
end

function PlotObservedVsModelMatrices(metabs, RModel, RObserved, RTrue, scenarioName, outputFile)

    fig = figure("Name", "Observed vs model: " + scenarioName, "Color", "w", "Position", [100, 100, 1400, 450]);
    tiledlayout(1, 3, "Padding", "compact", "TileSpacing", "compact");

    ax = nexttile;
    PlotHeatmap(ax, RModel, metabs, "Model R", [-1, 1], "gray");

    ax = nexttile;
    PlotHeatmap(ax, RObserved, metabs, "Mean empirical R", [-1, 1], "gray");

    ax = nexttile;
    PlotHeatmap(ax, RObserved - RModel, metabs, "Mean empirical R - model R", [-1, 1], "gray");

    titleText = "Scenario: " + scenarioName + " | true changed pairs are in R_{true} - R_{model}";
    sgtitle(titleText, ...
        "Interpreter", "none", ...
        "Color", [0 0 0], ...
        "FontWeight", "bold");

    SaveFigure(fig, outputFile);
end

function PlotHeatmap(ax, M, labels, titleText, climVals, cmapType)

    imagesc(ax, M);
    axis(ax, "image");

    cb = colorbar(ax);              % changed: keep colorbar handle
    caxis(ax, climVals);

    labels = string(labels(:));

    xticks(ax, 1:numel(labels));
    xticklabels(ax, labels);
    xtickangle(ax, 45);

    yticks(ax, 1:numel(labels));
    yticklabels(ax, labels);

    % Explicit non-transparent / high-contrast axes styling
    set(ax, ...
        "TickLabelInterpreter", "none", ...
        "FontSize", 12, ...
        "FontWeight", "bold", ...
        "XColor", [0 0 0], ...
        "YColor", [0 0 0], ...
        "Color", [1 1 1]);

    % Explicit non-transparent / high-contrast title styling
    title(ax, titleText, ...
        "Interpreter", "none", ...
        "Color", [0 0 0], ...
        "FontWeight", "bold");

    % Explicit colorbar styling
    cb.Color = [0 0 0];
    cb.FontWeight = "bold";

    if strcmpi(cmapType, "diverging")
        colormap(ax, BlueWhiteRed(256));
    else
        colormap(ax, gray(256));
    end

    for i = 1:size(M, 1)
        for j = 1:size(M, 2)
            text(ax, j, i, sprintf("%.2f", M(i, j)), ...
                "HorizontalAlignment", "center", ...
                "FontSize", 10, ...
                "FontWeight", "bold", ...
                "Color", [0 0 0]);
        end
    end
end

function SaveFigure(fig, outputFile)

    [folderPath, ~, ~] = fileparts(outputFile);

    if ~isfolder(folderPath)
        mkdir(folderPath);
    end

    try
        exportgraphics(fig, outputFile, "Resolution", 200);
    catch
        saveas(fig, outputFile);
    end
end

function [capturedText, out] = CaptureFunctionCall(funcHandle)

    capturedText = evalc("out = funcHandle();");
end

function T = MatrixToSquareTable(M, labels)

    labels = string(labels(:));

    T = array2table(M, ...
        "VariableNames", cellstr(labels), ...
        "RowNames", cellstr(labels));
end

function R = CovToCorr(Sigma)

    sd = sqrt(diag(Sigma));
    R = Sigma ./ (sd * sd.');
    R(1:size(R, 1)+1:end) = 1;
end

function Rmean = FisherMeanCorrStack(Rstack)

    Rclamped = max(min(Rstack, 0.999999), -0.999999);
    Z = atanh(Rclamped);

    Zmean = mean(Z, 3, "omitnan");
    Rmean = tanh(Zmean);

    n = size(Rmean, 1);
    Rmean(1:n+1:end) = 1;
end

function baseMu = BuildBaseMeans(metabs)

    metabs = string(metabs(:));
    baseMu = nan(1, numel(metabs));

    for k = 1:numel(metabs)
        switch metabs(k)
            case "NAA"
                baseMu(k) = 10.0;
            case "Cr"
                baseMu(k) = 7.0;
            case "PCr"
                baseMu(k) = 5.0;
            case "Glu"
                baseMu(k) = 12.0;
            otherwise
                baseMu(k) = 6.0 + k;
        end
    end
end

function idx = FindMetabIndex(metabs, target)

    metabs = string(metabs(:));
    target = string(target);

    idx = find(metabs == target, 1);

    if isempty(idx)
        error("Metabolite not found in toy metabolite list: %s", target);
    end
end

function AssertPositiveDefinite(M, matrixName)

    [~, flag] = chol(M);

    if flag ~= 0
        error("%s is not positive definite.", matrixName);
    end
end

function Q = BuildToyCRLBSummary(metabs, nParts)

    metabs = string(metabs(:));
    n = numel(metabs);

    Q = table( ...
        metabs, ...
        10 * ones(n, 1), ...
        nParts * ones(n, 1), ...
        nParts * ones(n, 1), ...
        ones(n, 1), ...
        false(n, 1), ...
        'VariableNames', { ...
            'metabolite', ...
            'medianCRLB', ...
            'nCRLBUnder100', ...
            'nInstances', ...
            'fractionCRLBUnder100', ...
            'fails90PercentRule'});
end

function Q = BuildToyCRLBQuality(metabs, nPatients, nParts)

    metabs = string(metabs(:));
    n = numel(metabs);
    nInstances = nPatients * nParts;

    Q = table( ...
        metabs, ...
        nInstances * ones(n, 1), ...
        nInstances * ones(n, 1), ...
        ones(n, 1), ...
        false(n, 1), ...
        'VariableNames', { ...
            'metabolite', ...
            'nCRLBUnder100', ...
            'nInstances', ...
            'fractionCRLBUnder100', ...
            'fails90PercentRule'});
end

function value = ExtractNumeric(T, varName, defaultValue)

    if nargin < 3
        defaultValue = NaN;
    end

    varName = string(varName);

    if isempty(T) || height(T) == 0 || ~ismember(varName, string(T.Properties.VariableNames))
        value = defaultValue;
        return;
    end

    value = T{1, varName};

    if iscell(value)
        value = value{1};
    end

    value = double(value);
end

function value = ExtractLogical(T, varName, defaultValue)

    if nargin < 3
        defaultValue = false;
    end

    varName = string(varName);

    if isempty(T) || height(T) == 0 || ~ismember(varName, string(T.Properties.VariableNames))
        value = defaultValue;
        return;
    end

    value = T{1, varName};

    if iscell(value)
        value = value{1};
    end

    value = logical(value);
end

function WriteTextFile(filePath, txt)

    fid = fopen(filePath, "w");

    if fid < 0
        warning("Could not write text file: %s", filePath);
        return;
    end

    cleaner = onCleanup(@() fclose(fid));
    fprintf(fid, "%s", txt);
end
