cfg = ProjectConfig();
cfg.covariance.loadMode = "allSubfolders";
cfg.covariance.ignoreZeros = false;
covOutputs = MetabCovarianceByPatient(cfg);

cfg.degraaf.loadMode = "allSubfolders";
cfg.degraaf.division = 1;
cfg.degraaf.maskInvalidCRLB = false;
cfg.degraaf.invalidCRLBValue = 100;
deGraafOutputs = DeGraafAmplitudeCorrelationByPatient(cfg);

filterCfg = struct();
filterCfg.patientIDs = "all";
filterCfg.division = 1;
filterCfg.metabolites = "all";
filterCfg.useSumPreferredFilter = true;
filterCfg.sumMetabolites = ["GPC+PCh", "NAA+NAAG", "Cr+PCr", "Glu+Gln"];
filterCfg.useCRLBMajorityFilter = true;
filterCfg.crlbMajorityThreshold = 100;
filterCfg.ignoreZeros = true;
filterCfg.pairwiseMinValidParts = 10;
filterCfg.pairwiseMinPatients = 3;
filterCfg.temporalMinValidParts = 8;
filterCfg.temporalMinPatients = 3;
filterCfg.temporalUseGlobalMetabolites = false;
filterCfg.temporalCRLBThreshold = 100;
filterCfg.temporalRequiredGoodFraction = 0.01;
filterCfg.prepareTemporalCircularShift = false;
filterCfg.wishartMinValidParts = 30;
filterCfg.wishartViews.modeA = struct('metabolites', "all", 'minValidParts', 30);
filterCfg.wishartViews.modeB = struct( ...
    'metabolites', ["NAA", "Cr", "PCr", "Glu"], 'minValidParts', 30);
filterCfg.wishartViews.modeC = struct('metabolites', "all", 'minValidParts', 30);
[analysisData, filterReport] = ApplyAnalysisFilters( ...
    covOutputs, deGraafOutputs, filterCfg); %#ok<ASGLU>

statsCfg = struct('alpha', 0.05, 'useFDR', true, 'exportResults', false);
basisCfg = struct('makeFigures', true, 'exportResults', true);
basisPairDiagnostics = AnalyzePostFitBasisPairwise( ...
    covOutputs, deGraafOutputs, analysisData, statsCfg, cfg, basisCfg);
save(fullfile(cfg.paths.rootDir, "PairwiseBasisModel_Comparison", ...
    "BasisPairDiagnostics.mat"), 'basisPairDiagnostics', '-v7.3');
