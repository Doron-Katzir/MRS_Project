%% Example: pairwise empirical-vs-LCModel/de Graaf correlation test
%
% Put this after you already created:
%   covOutputs = MetabCovarianceByPatient(cfg);
%   deGraafOutputs = DeGraafAmplitudeCorrelationByPatient(cfg);

pairCfg = struct();

% Use all patients present in both outputs.
pairCfg.patientIDs = "all";

% Choose metabolites to test.
% Option 1: all available non-sum metabolites
pairCfg.metabolites = "all";

% Option 2: manually specify a clean panel
% pairCfg.metabolites = ["NAA", "Cr", "PCr", "Glu"];

% Exclude summed metabolites by default.
pairCfg.excludeSumMetabolites = true;
pairCfg.sumMetabolites = ["GPC+PCh", "NAA+NAAG", "Cr+PCr", "Glu+Gln"];

% Pairwise correlations are computed using valid repeated parts for that pair.
pairCfg.ignoreZeros = true;
pairCfg.minValidParts = 10;

% Minimum number of patients required to run a group-level t-test for a pair.
pairCfg.minPatientsForGroupTest = 3;

% Significance threshold and export.
pairCfg.alpha = 0.05;
pairCfg.exportResults = true;
pairCfg.outputDir = fullfile(pwd, "PairwiseEmpiricalVsModelResults");

pairOutputs = TestPairwiseEmpiricalVsModelCorrelation(covOutputs, deGraafOutputs, pairCfg);

%% Display main results

disp("Pairwise empirical-vs-LCModel/deGraaf summary:")
disp(pairOutputs.pairSummaryTable)

disp("Patient-level pairwise results:")
disp(pairOutputs.patientPairTable)

%% Example: show only FDR-significant pairs

sigPairs = pairOutputs.pairSummaryTable(pairOutputs.pairSummaryTable.rejectH0_FDR == true, :);

disp("FDR-significant empirical-vs-model pairs:")
disp(sigPairs)
