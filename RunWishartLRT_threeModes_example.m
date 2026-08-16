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
