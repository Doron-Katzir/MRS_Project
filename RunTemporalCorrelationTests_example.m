%% Example caller for TestTemporalMetaboliteCorrelations.m
%
% Run this after:
%   covOutputs = MetabCovarianceByPatient(cfg);
%   deGraafOutputs = DeGraafAmplitudeCorrelationByPatient(cfg);

%% Configure statistical tests

testCfg = struct();

% Patients to include:
%   "all"
%   ["P01", "P02", "P11"]
%   1:10
testCfg.patientIDs = "all";

% Metabolites for the full Fisher-z group test:
%   "all" tests all metabolites in covOutputs.metabList
%   or provide a selected list
testCfg.metabList = "all";
% Example subset:
% testCfg.metabList = ["Glu", "Gln", "Glu+Gln", "GABA", "NAA+NAAG", "Cr+PCr", "GPC+PCh"];

% Selected pairs for circular-shift permutation test.
% This test is more expensive, so use it on pairs you care about most.
testCfg.permutationPairs = [
    "Glu",     "GABA"
    "Glu+Gln", "GABA"
    "Glu",     "Gln"
    "Glu+Gln", "NAA+NAAG"
    "GPC+PCh", "Cr+PCr"
];

% Minimum number of valid Division_1 parts required inside one patient
testCfg.minValidParts = 8;

% Minimum number of patients required for group-level test
testCfg.minPatients = 3;

% CRLB reliability rule
testCfg.crlbThreshold = 100;
testCfg.requiredGoodFraction = 0.90;

% Which tests to run
testCfg.doFisherGroupTest = true;
testCfg.doCircularShiftPermutation = true;

% Circular-shift group permutation settings
testCfg.nGroupPermutations = 5000;
testCfg.rngSeed = 1;

% Export CSV tables
testCfg.exportTables = true;
testCfg.outputDir = fullfile(pwd, "TemporalCorrelationStats");

%% Run tests

testOutputs = TestTemporalMetaboliteCorrelations( ...
    covOutputs, ...
    deGraafOutputs, ...
    testCfg);

%% Inspect main outputs

% Main group-level statistical table:
testOutputs.groupFisherTable

% One row per patient per metabolite pair:
testOutputs.patientCorrelationTable

% Circular-shift temporal robustness test:
testOutputs.circularShiftTable

% CRLB reliability table used for pair status:
testOutputs.groupCRLBQualityTable

%% Useful quick views

% Top Fisher-z group-test results by FDR q-value:
sortrows(testOutputs.groupFisherTable, "qValue_FDR", "ascend")

% Pairs significant at FDR q < 0.05:
testOutputs.groupFisherTable(testOutputs.groupFisherTable.qValue_FDR < 0.05, :)

% Pairs that are significant and CRLB reliable:
idxReliableSignificant = ...
    testOutputs.groupFisherTable.qValue_FDR < 0.05 & ...
    testOutputs.groupFisherTable.CRLB_pair_status == "PASS";

testOutputs.groupFisherTable(idxReliableSignificant, :)

% Circular-shift test results sorted by q-value:
sortrows(testOutputs.circularShiftTable, "qValue_FDR", "ascend")
