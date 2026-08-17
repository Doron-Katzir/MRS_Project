function testOutputs = TestTemporalMetaboliteCorrelations(covOutputs, deGraafOutputs, testCfg)
% TestTemporalMetaboliteCorrelations
%
% Statistical tests for temporal metabolite correlations.
%
% This function implements two complementary tests:
%
%   TEST 1: Patient-level Fisher-z group test
%       For each metabolite pair:
%           1. Compute Pearson correlation within each patient across Division_1 parts.
%           2. Fisher-transform patient correlations: z = atanh(r).
%           3. Test across patients whether mean(z) differs from zero.
%           4. Correct p-values across metabolite pairs using Benjamini-Hochberg FDR.
%
%   TEST 2: Circular-shift permutation test
%       For selected metabolite pairs:
%           1. Compute zero-lag correlation within each patient.
%           2. Circularly shift one metabolite relative to the other.
%           3. Recompute correlations under shifted alignments.
%           4. Build a group-level null distribution by sampling one shift per patient.
%           5. Ask whether the observed zero-lag group correlation is stronger than expected.
%
% Inputs
% ------
% covOutputs:
%   Output from MetabCovarianceByPatient(cfg)
%
% deGraafOutputs:
%   Output from DeGraafAmplitudeCorrelationByPatient(cfg)
%   Use [] if unavailable.
%
% testCfg:
%   Struct with optional fields:
%       patientIDs                  "all", string list, or numeric indices
%       metabList                   "all" or string list
%       permutationPairs            N x 2 string array of selected pairs
%       minValidParts               minimum valid parts for within-patient correlation
%       minPatients                 minimum patients for group-level test
%       crlbThreshold               CRLB threshold, default 100
%       requiredGoodFraction        reliability threshold, default 0.90
%       doFisherGroupTest           true/false
%       doCircularShiftPermutation  true/false
%       nGroupPermutations          number of group-level permutation draws
%       rngSeed                     random seed
%       exportTables                true/false
%       outputDir                   output directory for CSV files
%
% Output
% ------
% testOutputs.groupFisherTable
% testOutputs.patientCorrelationTable
% testOutputs.circularShiftTable
% testOutputs.groupCRLBQualityTable
%
% Example
% -------
% testCfg = struct();
% testCfg.patientIDs = "all";
% testCfg.metabList = ["Glu", "Gln", "Glu+Gln", "GABA", "NAA+NAAG", "Cr+PCr"];
% testCfg.permutationPairs = ["Glu", "GABA"; "Glu+Gln", "GABA"];
% testOutputs = TestTemporalMetaboliteCorrelations(covOutputs, deGraafOutputs, testCfg);

if nargin < 2
    deGraafOutputs = [];
end

if nargin < 3 || isempty(testCfg)
    testCfg = struct();
end

testCfg = ApplyTestDefaults(testCfg);

if testCfg.exportTables && ~isfolder(testCfg.outputDir)
    mkdir(testCfg.outputDir);
end

rng(testCfg.rngSeed);

patientIDsAll = string(fieldnames(covOutputs.patientResultsByID));
patientIDs = ResolvePatientIDs(testCfg.patientIDs, patientIDsAll);
metabList = ResolveMetabList(testCfg.metabList, covOutputs);

testOutputs = struct();
testOutputs.config = testCfg;
testOutputs.patientIDs = patientIDs;
testOutputs.metabList = metabList;

fprintf('\nTemporal metabolite-correlation statistical tests\n');
fprintf('Patients selected: %d\n', numel(patientIDs));
fprintf('Metabolites selected: %d\n', numel(metabList));
fprintf('Minimum valid parts per patient: %d\n', testCfg.minValidParts);
fprintf('Minimum patients per group test: %d\n', testCfg.minPatients);

%% Build group-level CRLB reliability table

groupCRLBQualityTable = BuildGroupCRLBQualityTable( ...
    covOutputs, ...
    patientIDs, ...
    metabList, ...
    testCfg.crlbThreshold, ...
    testCfg.requiredGoodFraction);

testOutputs.groupCRLBQualityTable = groupCRLBQualityTable;

%% Test 1: patient-level Fisher-z group test

if testCfg.doFisherGroupTest

    fprintf('\nRunning patient-level Fisher-z group test...\n');

    [groupFisherTable, patientCorrelationTable] = RunFisherGroupTest( ...
        covOutputs, ...
        deGraafOutputs, ...
        patientIDs, ...
        metabList, ...
        groupCRLBQualityTable, ...
        testCfg);

    testOutputs.groupFisherTable = groupFisherTable;
    testOutputs.patientCorrelationTable = patientCorrelationTable;

    if testCfg.exportTables

        writetable(groupFisherTable, ...
            fullfile(testCfg.outputDir, "Group_Fisher_Z_Tests.csv"));

        writetable(patientCorrelationTable, ...
            fullfile(testCfg.outputDir, "Patient_Level_Correlations.csv"));
    end
end

%% Test 2: circular-shift permutation test

if testCfg.doCircularShiftPermutation

    fprintf('\nRunning circular-shift permutation test...\n');

    circularShiftTable = RunCircularShiftPermutationTest( ...
        covOutputs, ...
        patientIDs, ...
        testCfg.permutationPairs, ...
        testCfg);

    testOutputs.circularShiftTable = circularShiftTable;

    if testCfg.exportTables
        writetable(circularShiftTable, ...
            fullfile(testCfg.outputDir, "Circular_Shift_Permutation_Tests.csv"));
    end
end

if testCfg.exportTables
    writetable(groupCRLBQualityTable, ...
        fullfile(testCfg.outputDir, "CRLB_Reliability_Table.csv"));

    fprintf('\nSaved statistical-test tables to:\n%s\n', testCfg.outputDir);
end

fprintf('\nFinished temporal correlation statistical tests.\n');

end

%% ========================================================================
% Defaults
% ========================================================================

function testCfg = ApplyTestDefaults(testCfg)

testCfg = SetDefault(testCfg, "patientIDs", "all");
testCfg = SetDefault(testCfg, "metabList", "all");

testCfg = SetDefault(testCfg, "permutationPairs", [ ...
    "Glu", "GABA"; ...
    "Glu+Gln", "GABA"; ...
    "Glu", "Gln"; ...
    "Glu+Gln", "NAA+NAAG"; ...
    "GPC+PCh", "Cr+PCr"]);

testCfg = SetDefault(testCfg, "minValidParts", 8);
testCfg = SetDefault(testCfg, "minPatients", 3);

testCfg = SetDefault(testCfg, "crlbThreshold", 100);
testCfg = SetDefault(testCfg, "requiredGoodFraction", 0.90);

testCfg = SetDefault(testCfg, "doFisherGroupTest", true);
testCfg = SetDefault(testCfg, "doCircularShiftPermutation", true);

testCfg = SetDefault(testCfg, "nGroupPermutations", 5000);
testCfg = SetDefault(testCfg, "rngSeed", 1);

testCfg = SetDefault(testCfg, "exportTables", true);
testCfg = SetDefault(testCfg, "outputDir", fullfile(pwd, "TemporalCorrelationStats"));

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
    error("Unsupported testCfg.patientIDs format.");
end

missing = patientIDs(~ismember(patientIDs, patientIDsAll));

if ~isempty(missing)
    error("Requested patientIDs were not found: %s", strjoin(missing, ", "));
end

end

function metabList = ResolveMetabList(metabOption, covOutputs)

if ischar(metabOption) || isstring(metabOption)

    metabOption = string(metabOption);

    if isscalar(metabOption) && metabOption == "all"
        metabList = string(covOutputs.metabList(:));
    else
        metabList = metabOption(:);
    end

else
    metabList = string(metabOption(:));
end

end

%% ========================================================================
% Test 1: patient-level Fisher-z group test
% ========================================================================

function [groupTable, patientTable] = RunFisherGroupTest(covOutputs, deGraafOutputs, patientIDs, metabList, crlbQualityTable, testCfg)

nMetabs = numel(metabList);

pairRows = table();
patientRows = table();

for a = 1:nMetabs-1

    for b = a+1:nMetabs

        metabA = metabList(a);
        metabB = metabList(b);

        patientR = nan(numel(patientIDs), 1);
        patientZ = nan(numel(patientIDs), 1);
        patientN = nan(numel(patientIDs), 1);
        patientUsed = false(numel(patientIDs), 1);

        for pIdx = 1:numel(patientIDs)

            patientID = patientIDs(pIdx);
            T = covOutputs.patientResultsByID.(char(patientID)).partTable;

            [rVal, nValid] = ComputePatientCorrelation(T, metabA, metabB, testCfg.minValidParts);

            patientR(pIdx) = rVal;
            patientN(pIdx) = nValid;

            if ~isnan(rVal)
                patientZ(pIdx) = FisherZ(rVal);
                patientUsed(pIdx) = true;
            end

            patientRows = [patientRows; table( ... %#ok<AGROW>
                patientID, ...
                metabA, ...
                metabB, ...
                rVal, ...
                patientZ(pIdx), ...
                nValid, ...
                patientUsed(pIdx), ...
                'VariableNames', { ...
                'patientID', ...
                'metaboliteA', ...
                'metaboliteB', ...
                'r', ...
                'fisherZ', ...
                'nValidParts', ...
                'usedInGroupTest'})];
        end

        valid = ~isnan(patientZ);
        zVals = patientZ(valid);
        rVals = patientR(valid);

        nPatientsUsed = sum(valid);

        if nPatientsUsed >= testCfg.minPatients

            meanZ = mean(zVals);
            sdZ = std(zVals, 0);
            seZ = sdZ / sqrt(nPatientsUsed);
            groupMeanR = tanh(meanZ);
            medianR = median(rVals);

            if sdZ == 0
                if meanZ == 0
                    tStat = 0;
                    pValue = 1;
                else
                    tStat = sign(meanZ) * Inf;
                    pValue = 0;
                end
            else
                tStat = meanZ / seZ;
                df = nPatientsUsed - 1;
                pValue = 2 * (1 - StudentTCDF(abs(tStat), df));
            end

        else

            meanZ = NaN;
            sdZ = NaN;
            seZ = NaN;
            groupMeanR = NaN;
            medianR = NaN;
            tStat = NaN;
            df = NaN;
            pValue = NaN;
        end

        nPositive = sum(rVals > 0);
        nNegative = sum(rVals < 0);
        nZero = sum(rVals == 0);

        if nPatientsUsed > 0
            signConsistency = max(nPositive, nNegative) / nPatientsUsed;
        else
            signConsistency = NaN;
        end

        [fracA, fracB, pairCRLBStatus] = GetGroupCRLBStatus( ...
            metabA, metabB, crlbQualityTable, testCfg.requiredGoodFraction);

        [signedLCModelCorr, absLCModelCorr] = LookupDeGraafCorrelation( ...
            deGraafOutputs, metabA, metabB);

        [absEmpiricalCorrFromMatrix] = LookupEmpiricalAbsCorrelation( ...
            covOutputs, metabA, metabB);

        diff_LCM_minus_empirical = absLCModelCorr - absEmpiricalCorrFromMatrix;

        if ~isnan(absLCModelCorr) && absLCModelCorr ~= 0
            empirical_over_LCModel = absEmpiricalCorrFromMatrix / absLCModelCorr;
            relativeDiff_LCM_minus_empirical_percent = ...
                100 * (absLCModelCorr - absEmpiricalCorrFromMatrix) / absLCModelCorr;
        else
            empirical_over_LCModel = NaN;
            relativeDiff_LCM_minus_empirical_percent = NaN;
        end

        pairRows = [pairRows; table( ... %#ok<AGROW>
            metabA, ...
            metabB, ...
            nPatientsUsed, ...
            groupMeanR, ...
            medianR, ...
            meanZ, ...
            sdZ, ...
            seZ, ...
            tStat, ...
            df, ...
            pValue, ...
            nPositive, ...
            nNegative, ...
            nZero, ...
            signConsistency, ...
            pairCRLBStatus, ...
            fracA, ...
            fracB, ...
            signedLCModelCorr, ...
            absLCModelCorr, ...
            absEmpiricalCorrFromMatrix, ...
            diff_LCM_minus_empirical, ...
            empirical_over_LCModel, ...
            relativeDiff_LCM_minus_empirical_percent, ...
            'VariableNames', { ...
            'metaboliteA', ...
            'metaboliteB', ...
            'nPatientsUsed', ...
            'groupMeanR', ...
            'medianPatientR', ...
            'meanFisherZ', ...
            'sdFisherZ', ...
            'seFisherZ', ...
            'tStat', ...
            'df', ...
            'pValue', ...
            'nPositivePatients', ...
            'nNegativePatients', ...
            'nZeroPatients', ...
            'signConsistency', ...
            'CRLB_pair_status', ...
            'metaboliteA_fractionCRLBUnder100', ...
            'metaboliteB_fractionCRLBUnder100', ...
            'signedLCModelAmplitudeCorr', ...
            'absLCModelAmplitudeCorr', ...
            'absEmpiricalCorrFromMatrix', ...
            'diff_absLCModel_minus_absEmpirical', ...
            'empirical_over_LCModel', ...
            'relativeDiff_LCM_minus_empirical_percent'})];
    end
end

groupTable = pairRows;
groupTable.qValue_FDR = BenjaminiHochbergFDR(groupTable.pValue);

groupTable = sortrows(groupTable, "pValue", "ascend");
patientTable = patientRows;

end

function [rVal, nValid] = ComputePatientCorrelation(T, metabA, metabB, minValidParts)

rVal = NaN;
nValid = 0;

T = SortPartTable(T);

colA = matlab.lang.makeValidName(char(metabA));
colB = matlab.lang.makeValidName(char(metabB));

if ~ismember(colA, string(T.Properties.VariableNames)) || ...
   ~ismember(colB, string(T.Properties.VariableNames))
    return;
end

x = T.(colA);
y = T.(colB);

valid = ~isnan(x) & ~isnan(y);
nValid = sum(valid);

if nValid < minValidParts
    return;
end

rVal = corr(x(valid), y(valid));

end

function z = FisherZ(r)

r = max(min(r, 0.999999), -0.999999);
z = atanh(r);

end

%% ========================================================================
% Test 2: circular-shift permutation test
% ========================================================================

function circularShiftTable = RunCircularShiftPermutationTest(covOutputs, patientIDs, permutationPairs, testCfg)

if isempty(permutationPairs)
    circularShiftTable = table();
    return;
end

permutationPairs = string(permutationPairs);

nPairs = size(permutationPairs, 1);
resultRows = table();

for pairIdx = 1:nPairs

    metabA = permutationPairs(pairIdx, 1);
    metabB = permutationPairs(pairIdx, 2);

    observedZ = nan(numel(patientIDs), 1);
    observedR = nan(numel(patientIDs), 1);
    perPatientP = nan(numel(patientIDs), 1);
    perPatientN = nan(numel(patientIDs), 1);
    shiftZDistributions = cell(numel(patientIDs), 1);

    for pIdx = 1:numel(patientIDs)

        patientID = patientIDs(pIdx);
        T = covOutputs.patientResultsByID.(char(patientID)).partTable;

        [rObs, nObs, rShiftVals] = ComputeCircularShiftCorrelationsForPatient( ...
            T, metabA, metabB, testCfg.minValidParts);

        observedR(pIdx) = rObs;
        perPatientN(pIdx) = nObs;

        if ~isnan(rObs)
            observedZ(pIdx) = FisherZ(rObs);
        end

        if ~isempty(rShiftVals) && ~isnan(rObs)
            perPatientP(pIdx) = ...
                (1 + sum(abs(rShiftVals) >= abs(rObs))) / (1 + numel(rShiftVals));
            shiftZDistributions{pIdx} = FisherZ(rShiftVals(:));
        end
    end

    validPatients = ~isnan(observedZ) & CellHasValues(shiftZDistributions);
    nPatientsUsed = sum(validPatients);

    if nPatientsUsed >= testCfg.minPatients

        observedGroupMeanZ = mean(observedZ(validPatients));
        observedGroupR = tanh(observedGroupMeanZ);

        permMeanZ = nan(testCfg.nGroupPermutations, 1);
        validIdx = find(validPatients);

        for permIdx = 1:testCfg.nGroupPermutations

            curZ = nan(numel(validIdx), 1);

            for j = 1:numel(validIdx)

                pIdx = validIdx(j);
                zPool = shiftZDistributions{pIdx};
                randomIndex = randi(numel(zPool));
                curZ(j) = zPool(randomIndex);
            end

            permMeanZ(permIdx) = mean(curZ);
        end

        groupPermutationP = ...
            (1 + sum(abs(permMeanZ) >= abs(observedGroupMeanZ))) / ...
            (1 + numel(permMeanZ));

        meanPatientR = mean(observedR(validPatients));
        medianPatientR = median(observedR(validPatients));
        medianPatientPermutationP = median(perPatientP(validPatients), 'omitnan');

        nPositive = sum(observedR(validPatients) > 0);
        nNegative = sum(observedR(validPatients) < 0);
        signConsistency = max(nPositive, nNegative) / nPatientsUsed;

    else

        observedGroupMeanZ = NaN;
        observedGroupR = NaN;
        groupPermutationP = NaN;
        meanPatientR = NaN;
        medianPatientR = NaN;
        medianPatientPermutationP = NaN;
        nPositive = NaN;
        nNegative = NaN;
        signConsistency = NaN;
    end

    resultRows = [resultRows; table( ... %#ok<AGROW>
        metabA, ...
        metabB, ...
        nPatientsUsed, ...
        observedGroupR, ...
        observedGroupMeanZ, ...
        meanPatientR, ...
        medianPatientR, ...
        groupPermutationP, ...
        medianPatientPermutationP, ...
        nPositive, ...
        nNegative, ...
        signConsistency, ...
        'VariableNames', { ...
        'metaboliteA', ...
        'metaboliteB', ...
        'nPatientsUsed', ...
        'observedGroupR', ...
        'observedGroupMeanFisherZ', ...
        'meanPatientR', ...
        'medianPatientR', ...
        'groupCircularShiftPValue', ...
        'medianPatientCircularShiftPValue', ...
        'nPositivePatients', ...
        'nNegativePatients', ...
        'signConsistency'})];
end

circularShiftTable = resultRows;
circularShiftTable.qValue_FDR = BenjaminiHochbergFDR(circularShiftTable.groupCircularShiftPValue);
circularShiftTable = sortrows(circularShiftTable, "groupCircularShiftPValue", "ascend");

end

function [rObs, nObs, rShiftVals] = ComputeCircularShiftCorrelationsForPatient(T, metabA, metabB, minValidParts)

rObs = NaN;
nObs = 0;
rShiftVals = [];

T = SortPartTable(T);

colA = matlab.lang.makeValidName(char(metabA));
colB = matlab.lang.makeValidName(char(metabB));

if ~ismember(colA, string(T.Properties.VariableNames)) || ...
   ~ismember(colB, string(T.Properties.VariableNames))
    return;
end

x = T.(colA);
y = T.(colB);

valid = ~isnan(x) & ~isnan(y);
nObs = sum(valid);

if nObs < minValidParts
    return;
end

rObs = corr(x(valid), y(valid));

nTime = numel(x);
rShiftVals = nan(nTime - 1, 1);

for shiftVal = 1:nTime-1

    yShift = circshift(y, shiftVal);
    validShift = ~isnan(x) & ~isnan(yShift);

    if sum(validShift) >= minValidParts
        rShiftVals(shiftVal) = corr(x(validShift), yShift(validShift));
    end
end

rShiftVals = rShiftVals(~isnan(rShiftVals));

end

function mask = CellHasValues(C)

mask = false(numel(C), 1);

for i = 1:numel(C)
    mask(i) = ~isempty(C{i}) && any(~isnan(C{i}));
end

end

%% ========================================================================
% CRLB reliability table
% ========================================================================

function crlbQualityTable = BuildGroupCRLBQualityTable(covOutputs, patientIDs, metabList, crlbThreshold, requiredGoodFraction)

metabList = string(metabList(:));
nMetabs = numel(metabList);

goodCRLBCount = zeros(nMetabs, 1);
totalInstanceCount = zeros(nMetabs, 1);

for pIdx = 1:numel(patientIDs)

    patientID = patientIDs(pIdx);

    partTable = covOutputs.patientResultsByID.(char(patientID)).partTable;
    coordTable = covOutputs.patientResultsByID.(char(patientID)).coordTable;

    patientTable = BuildPatientCRLBQualityTable( ...
        covOutputs, ...
        metabList, ...
        coordTable, ...
        partTable, ...
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

function patientCRLBQualityTable = BuildPatientCRLBQualityTable(covOutputs, metabList, coordTable, partTable, crlbThreshold, requiredGoodFraction)

metabList = string(metabList(:));
nMetabs = numel(metabList);

if isfield(covOutputs, "settings") && isfield(covOutputs.settings, "division")
    divisionUsed = covOutputs.settings.division;
else
    divisionUsed = 1;
end

partsUsed = partTable.part;

coordTable.name = string(coordTable.name);
coordTable.filename = string(coordTable.filename);

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

error("Could not find a CRLB / %%SD column. Available columns are:\n%s", ...
    strjoin(varNames, ", "));

end

function [fracA, fracB, status] = GetGroupCRLBStatus(metabA, metabB, crlbQualityTable, requiredGoodFraction)

crlbQualityTable.metabolite = string(crlbQualityTable.metabolite);

idxA = find(crlbQualityTable.metabolite == string(metabA), 1);
idxB = find(crlbQualityTable.metabolite == string(metabB), 1);

if isempty(idxA)
    fracA = NaN;
else
    fracA = crlbQualityTable.fractionCRLBUnder100(idxA);
end

if isempty(idxB)
    fracB = NaN;
else
    fracB = crlbQualityTable.fractionCRLBUnder100(idxB);
end

if fracA >= requiredGoodFraction && fracB >= requiredGoodFraction
    status = "PASS";
else
    status = "FAIL";
end

end

%% ========================================================================
% LCModel / De Graaf and empirical lookup helpers
% ========================================================================

function [signedCorr, absCorr] = LookupDeGraafCorrelation(deGraafOutputs, metabA, metabB)

signedCorr = NaN;
absCorr = NaN;

if isempty(deGraafOutputs) || ~isstruct(deGraafOutputs) || ~isfield(deGraafOutputs, "group")
    return;
end

if isfield(deGraafOutputs.group, "meanAmplitudeCorrTable")
    signedCorr = LookupMatrixEntry(deGraafOutputs.group.meanAmplitudeCorrTable, metabA, metabB);
end

if isfield(deGraafOutputs.group, "meanAbsAmplitudeCorrTable")
    absCorr = LookupMatrixEntry(deGraafOutputs.group.meanAbsAmplitudeCorrTable, metabA, metabB);
elseif ~isnan(signedCorr)
    absCorr = abs(signedCorr);
end

end

function absCorr = LookupEmpiricalAbsCorrelation(covOutputs, metabA, metabB)

absCorr = NaN;

if isfield(covOutputs, "group") && isfield(covOutputs.group, "meanAbsCorrTable")
    absCorr = LookupMatrixEntry(covOutputs.group.meanAbsCorrTable, metabA, metabB);
end

end

function value = LookupMatrixEntry(T, metabA, metabB)

value = NaN;

if isempty(T)
    return;
end

labels = string(T.Properties.RowNames);

idxA = find(labels == string(metabA), 1);
idxB = find(labels == string(metabB), 1);

if isempty(idxA) || isempty(idxB)
    return;
end

M = table2array(T);
value = M(idxA, idxB);

end

%% ========================================================================
% Statistics helpers
% ========================================================================

function q = BenjaminiHochbergFDR(p)

p = p(:);
q = nan(size(p));

valid = ~isnan(p);
pValid = p(valid);

if isempty(pValid)
    return;
end

[pSorted, sortIdx] = sort(pValid, 'ascend');
m = numel(pSorted);

qSorted = pSorted .* m ./ (1:m)';

for i = m-1:-1:1
    qSorted(i) = min(qSorted(i), qSorted(i+1));
end

qSorted(qSorted > 1) = 1;

validIdx = find(valid);
q(validIdx(sortIdx)) = qSorted;

end

function F = StudentTCDF(x, v)
% StudentTCDF
% Small helper to avoid dependency on tcdf.
% Computes CDF of Student's t distribution for scalar x and degrees of freedom v.

if isnan(x) || isnan(v) || v <= 0
    F = NaN;
    return;
end

if isinf(x)
    if x > 0
        F = 1;
    else
        F = 0;
    end
    return;
end

if x == 0
    F = 0.5;
    return;
end

z = v / (v + x^2);
I = betainc(z, v/2, 0.5);

if x > 0
    F = 1 - 0.5 * I;
else
    F = 0.5 * I;
end

end

%% ========================================================================
% Small utility
% ========================================================================

function T = SortPartTable(T)

if ismember("part", string(T.Properties.VariableNames))
    [~, order] = sort(T.part);
    T = T(order, :);
end

end
