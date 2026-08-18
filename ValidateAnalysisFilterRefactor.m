function results = ValidateAnalysisFilterRefactor()
%ValidateAnalysisFilterRefactor Compare centralized outputs with saved legacy CSVs.

projectDir = fileparts(mfilename('fullpath'));
rootDir = fileparts(projectDir);
cfg = ProjectConfig();
cfg.paths.rootDir = rootDir;
cfg.paths.coordDir = fullfile(rootDir, "LCMFit");
cfg.load.mode = "allSubfolders";
cfg.load.coordDir = cfg.paths.coordDir;
cfg.covariance.loadMode = "allSubfolders";
cfg.covariance.ignoreZeros = false;
cfg.degraaf.loadMode = "allSubfolders";
cfg.degraaf.division = 1;
cfg.degraaf.maskInvalidCRLB = false;
cfg.degraaf.invalidCRLBValue = 100;

covOutputs = MetabCovarianceByPatient(cfg);
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
    covOutputs, deGraafOutputs, filterCfg);

pairCfg = struct('alpha', 0.05, 'exportResults', false);
pairOutputs = TestPairwiseEmpiricalVsModelCorrelation(analysisData, pairCfg);

summaryFile = fullfile(rootDir, "PairwiseEmpiricalVsModelResults", ...
    "Pairwise_Empirical_vs_Model_Summary.csv");
patientFile = fullfile(rootDir, "PairwiseEmpiricalVsModelResults", ...
    "Pairwise_Empirical_vs_Model_PatientLevel.csv");
legacySummary = readtable(summaryFile, 'TextType', 'string');
legacyPatient = readtable(patientFile, 'TextType', 'string');

results = struct();
results.filterReport = filterReport;
results.pairSummary = ComparePairSummary(pairOutputs.pairSummaryTable, legacySummary);
results.patientLevel = ComparePatientTable(pairOutputs.patientPairTable, legacyPatient);

baseLrtCfg = struct('minMetabolites', 2, 'alpha', 0.05, ...
    'applyNumericalRidge', true, 'ridgeScale', 1e-8, 'exportResults', false);
lrtCfg = baseLrtCfg;
lrtCfg.metaboliteSelectionMode = "perPatientLargestValid";
lrtCfg.filterView = "modeA";
modeA = TestWishartCovarianceLRT(analysisData, lrtCfg);
lrtCfg = baseLrtCfg;
lrtCfg.metaboliteSelectionMode = "fixed";
lrtCfg.filterView = "modeB";
modeB = TestWishartCovarianceLRT(analysisData, lrtCfg);
lrtCfg = baseLrtCfg;
lrtCfg.metaboliteSelectionMode = "largestCommon";
lrtCfg.filterView = "modeC";
lrtCfg.runOnlyCommonPanelPatients = true;
modeC = TestWishartCovarianceLRT(analysisData, lrtCfg);
results.wishartModeA = CompareWishart(modeA, fullfile(rootDir, ...
    "WishartLRTResults_ModeA_PerPatientLargestValid"));
results.wishartModeB = CompareWishart(modeB, fullfile(rootDir, ...
    "WishartLRTResults_ModeB_PredefinedPanel"));
results.wishartModeC = CompareWishart(modeC, fullfile(rootDir, ...
    "WishartLRTResults_ModeC_LargestCommon"));

temporalCfg = struct();
temporalCfg.doFisherGroupTest = true;
temporalCfg.doCircularShiftPermutation = false;
temporalCfg.exportTables = false;
temporalCfg.rngSeed = 1;
temporalOutputs = TestTemporalMetaboliteCorrelations(analysisData, temporalCfg);
oldTemporalGroup = readtable(fullfile(rootDir, "TemporalCorrelationStats", ...
    "Group_Fisher_Z_Tests.csv"), 'TextType', 'string');
oldTemporalPatient = readtable(fullfile(rootDir, "TemporalCorrelationStats", ...
    "Patient_Level_Correlations.csv"), 'TextType', 'string');
results.temporalGroup = CompareByPair(temporalOutputs.groupFisherTable, ...
    oldTemporalGroup, ["nPatientsUsed","groupMeanR","meanFisherZ", ...
    "tStat","pValue","qValue_FDR"]);
results.temporalPatient = CompareByPatientPair(temporalOutputs.patientCorrelationTable, ...
    oldTemporalPatient, ["r","fisherZ","nValidParts"]);

fprintf('\nCentral filtering validation against saved pre-refactor CSVs\n');
disp(results.pairSummary)
disp(results.patientLevel)
disp(results.wishartModeA)
disp(results.wishartModeB)
disp(results.wishartModeC)
disp(results.wishartModeA.patient)
disp(results.wishartModeA.group)
disp(results.wishartModeB.patient)
disp(results.wishartModeB.group)
disp(results.wishartModeC.patient)
disp(results.wishartModeC.group)
disp(results.temporalGroup)
disp(results.temporalPatient)
fprintf('Patients: %d; starting metabolites: %d; final metabolites: %d; pairs: %d\n', ...
    filterReport.commonPatientCount, filterReport.startingMetaboliteCount, ...
    filterReport.finalMetaboliteCount, height(pairOutputs.pairSummaryTable));
save(fullfile(rootDir, "GeneratedCache", "FilterRefactorValidationResults.mat"), ...
    'results');
end

function comparison = CompareWishart(newOutputs, legacyDir)
oldPatient = readtable(fullfile(legacyDir, "Wishart_LRT_PerPatient.csv"), ...
    'TextType', 'string');
oldGroup = readtable(fullfile(legacyDir, "Wishart_LRT_Group.csv"), ...
    'TextType', 'string');
patientComparison = CompareByKey(newOutputs.patientSummaryTable, oldPatient, ...
    "patientID", ["nValidParts","nMetabolites","df","T","pValue", ...
    "ridgeEmpirical","ridgeModel"]);
groupVariables = ["nPatientsUsed","T_group","df_group","pValue_group"];
groupDiff = MaxDifferences(newOutputs.groupTable, oldGroup, groupVariables);
comparison = struct();
comparison.patient = patientComparison;
comparison.group = table(groupVariables(:), groupDiff, ...
    'VariableNames', {'variable','maxAbsoluteDifference'});
end

function comparison = CompareByPair(newTable, oldTable, variables)
comparison = CompareByKey(newTable, oldTable, ["metaboliteA","metaboliteB"], variables);
end

function comparison = CompareByPatientPair(newTable, oldTable, variables)
comparison = CompareByKey(newTable, oldTable, ...
    ["patientID","metaboliteA","metaboliteB"], variables);
end

function comparison = CompareByKey(newTable, oldTable, keyVariables, variables)
newKey = repmat("", height(newTable), 1);
oldKey = repmat("", height(oldTable), 1);
for k = 1:numel(keyVariables)
    newKey = newKey + "|" + string(newTable.(keyVariables(k)));
    oldKey = oldKey + "|" + string(oldTable.(keyVariables(k)));
end
[matched, oldIdx] = ismember(newKey, oldKey);
maxDiff = nan(numel(variables), 1);
for k = 1:numel(variables)
    if all(matched) && ismember(variables(k), string(oldTable.Properties.VariableNames))
        x = double(newTable.(variables(k)));
        y = double(oldTable.(variables(k))(oldIdx));
        finite = isfinite(x) & isfinite(y);
        if any(finite), maxDiff(k) = max(abs(x(finite)-y(finite))); else, maxDiff(k) = 0; end
    end
end
comparison = table(variables(:), maxDiff, ...
    repmat(all(matched) && height(newTable) == height(oldTable), numel(variables), 1), ...
    'VariableNames', {'variable','maxAbsoluteDifference','allRowsMatched'});
end

function maxDiff = MaxDifferences(newTable, oldTable, variables)
maxDiff = nan(numel(variables), 1);
for k = 1:numel(variables)
    x = double(newTable.(variables(k)));
    y = double(oldTable.(variables(k)));
    finite = isfinite(x) & isfinite(y);
    if any(finite), maxDiff(k) = max(abs(x(finite)-y(finite))); else, maxDiff(k) = 0; end
end
end

function comparison = ComparePairSummary(newTable, oldTable)
newKey = string(newTable.metaboliteA) + "|" + string(newTable.metaboliteB);
oldKey = string(oldTable.metaboliteA) + "|" + string(oldTable.metaboliteB);
[matched, oldIdx] = ismember(newKey, oldKey);

variables = ["nPatientsUsed", "meanDeltaZ_empMinusModel", "tValue", ...
    "pValue", "qValue_FDR", "rejectH0", "rejectH0_FDR", ...
    "meanAbsEmpiricalR", "meanAbsModelR", ...
    "meanDeltaAbsR_empMinusModel"];
maxDiff = nan(numel(variables), 1);
for k = 1:numel(variables)
    if all(matched) && ismember(variables(k), string(oldTable.Properties.VariableNames))
        x = double(newTable.(variables(k)));
        y = double(oldTable.(variables(k))(oldIdx));
        finite = isfinite(x) & isfinite(y);
        if any(finite), maxDiff(k) = max(abs(x(finite) - y(finite))); else, maxDiff(k) = 0; end
    end
end
comparison = table(variables(:), maxDiff, ...
    'VariableNames', {'variable','maxAbsoluteDifference'});
comparison.allPairsMatched = repmat(all(matched) && height(newTable) == height(oldTable), ...
    height(comparison), 1);
end

function comparison = ComparePatientTable(newTable, oldTable)
newKey = string(newTable.patientID) + "|" + string(newTable.metaboliteA) + ...
    "|" + string(newTable.metaboliteB);
oldKey = string(oldTable.patientID) + "|" + string(oldTable.metaboliteA) + ...
    "|" + string(oldTable.metaboliteB);
[matched, oldIdx] = ismember(newKey, oldKey);
variables = ["nValidParts", "rEmpirical", "rModel", "deltaZ_empMinusModel"];
maxDiff = nan(numel(variables), 1);
for k = 1:numel(variables)
    if all(matched)
        x = double(newTable.(variables(k)));
        y = double(oldTable.(variables(k))(oldIdx));
        finite = isfinite(x) & isfinite(y);
        if any(finite), maxDiff(k) = max(abs(x(finite) - y(finite))); else, maxDiff(k) = 0; end
    end
end
comparison = table(variables(:), maxDiff, ...
    'VariableNames', {'variable','maxAbsoluteDifference'});
comparison.allRowsMatched = repmat(all(matched) && height(newTable) == height(oldTable), ...
    height(comparison), 1);
end
