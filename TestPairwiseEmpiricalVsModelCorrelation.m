function pairOutputs = TestPairwiseEmpiricalVsModelCorrelation(covOutputs, deGraafOutputs, pairCfg)
% TestPairwiseEmpiricalVsModelCorrelation
%
% Pairwise test against the LCModel/de Graaf model correlation.
%
% For every metabolite pair and every patient:
%   r_empirical,p = corr(metabolite A across Division_1 parts,
%                        metabolite B across Division_1 parts)
%   r_model,p     = LCModel/de Graaf amplitude correlation for that pair
%
% Fisher transform:
%   z_empirical,p = atanh(r_empirical,p)
%   z_model,p     = atanh(r_model,p)
%
% Difference:
%   d_p = z_empirical,p - z_model,p
%
% Group-level null hypothesis for each metabolite pair:
%   H0: mean(d_p) = 0
%
% The test is a one-sample t-test across patients, keeping patient as the
% independent unit.
%
% Expected inputs:
%   covOutputs.patientResultsByID.<patientID>.partTable
%   deGraafOutputs.patientResultsByID.<patientID>.meanAmplitudeCorrTable
% or:
%   deGraafOutputs.patientResultsByID.<patientID>.meanAmplitudeCovTable
%
% Example:
%   pairCfg = struct();
%   pairCfg.patientIDs = "all";
%   pairCfg.metabolites = ["NAA", "Cr", "PCr", "Glu"];
%   pairCfg.minValidParts = 10;
%   pairOutputs = TestPairwiseEmpiricalVsModelCorrelation(covOutputs, deGraafOutputs, pairCfg);
%   disp(pairOutputs.pairSummaryTable)

if nargin < 3 || isempty(pairCfg)
    pairCfg = struct();
end

opts = ParsePairwiseOptions(pairCfg);
patientIDs = ResolvePatientIDs(covOutputs, deGraafOutputs, opts.patientIDs);
[metabs, crlbMajorityTable] = ResolveMetabolites(covOutputs, deGraafOutputs, patientIDs, opts);

fprintf('\nRunning pairwise empirical-vs-LCModel/deGraaf correlation tests...\n');
fprintf('Number of patients available: %d\n', numel(patientIDs));
fprintf('Number of metabolites requested: %d\n', numel(metabs));

pairSummaryTable = table();
patientPairTable = table();

for aIdx = 1:numel(metabs)
    for bIdx = aIdx+1:numel(metabs)

        metabA = string(metabs(aIdx));
        metabB = string(metabs(bIdx));

        thisPatientRows = table();

        for pIdx = 1:numel(patientIDs)

            patientID = string(patientIDs(pIdx));
            patientField = matlab.lang.makeValidName(char(patientID));

            row = MakeEmptyPatientPairRow(patientID, metabA, metabB);

            try
                covPatient = covOutputs.patientResultsByID.(patientField);
                modelPatient = deGraafOutputs.patientResultsByID.(patientField);

                if ~isfield(covPatient, 'partTable')
                    error('Missing partTable.');
                end

                partTable = covPatient.partTable;

                x = GetTableColumnByMetab(partTable, metabA);
                y = GetTableColumnByMetab(partTable, metabB);

                if opts.ignoreZeros
                    x(x == 0) = NaN;
                    y(y == 0) = NaN;
                end

                valid = isfinite(x) & isfinite(y);
                xv = x(valid);
                yv = y(valid);

                row.nValidParts = numel(xv);

                if row.nValidParts < opts.minValidParts
                    error('Only %d valid repeated parts; minimum requested is %d.', row.nValidParts, opts.minValidParts);
                end

                rEmp = PearsonR(xv, yv);
                rModel = GetModelCorrelation(modelPatient, metabA, metabB);

                if ~isfinite(rEmp)
                    error('Empirical correlation is non-finite.');
                end
                if ~isfinite(rModel)
                    error('Model correlation is non-finite.');
                end

                zEmp = atanh(ClampCorrelation(rEmp));
                zModel = atanh(ClampCorrelation(rModel));
                deltaZ = zEmp - zModel;

                row.status = "ok";
                row.rEmpirical = rEmp;
                row.rModel = rModel;
                row.zEmpirical = zEmp;
                row.zModel = zModel;
                row.deltaZ_empMinusModel = deltaZ;
                row.deltaR_empMinusModel = rEmp - rModel;
                row.errorMessage = "";

            catch ME
                row.status = "failed";
                row.errorMessage = string(ME.message);
            end

            thisPatientRows = [thisPatientRows; row]; %#ok<AGROW>
        end

        patientPairTable = [patientPairTable; thisPatientRows]; %#ok<AGROW>

        okRows = strcmp(thisPatientRows.status, "ok") & isfinite(thisPatientRows.deltaZ_empMinusModel);
        nPatientsUsed = sum(okRows);

        summaryRow = MakeEmptyPairSummaryRow(metabA, metabB);
        summaryRow.nPatientsUsed = nPatientsUsed;

        if nPatientsUsed >= opts.minPatientsForGroupTest

            d = thisPatientRows.deltaZ_empMinusModel(okRows);
            zEmp = thisPatientRows.zEmpirical(okRows);
            zModel = thisPatientRows.zModel(okRows);
            rEmp = thisPatientRows.rEmpirical(okRows);
            rModel = thisPatientRows.rModel(okRows);
            nValidParts = thisPatientRows.nValidParts(okRows);

            [tValue, pValue, df] = OneSampleTTestAgainstZero(d);

            summaryRow.status = "ok";
            summaryRow.df = df;
            summaryRow.meanDeltaZ_empMinusModel = mean(d);
            summaryRow.sdDeltaZ_empMinusModel = std(d, 0);
            summaryRow.tValue = tValue;
            summaryRow.pValue = pValue;
            summaryRow.rejectH0 = pValue < opts.alpha;
            summaryRow.groupEmpiricalR = tanh(mean(zEmp));
            summaryRow.groupModelR = tanh(mean(zModel));
            summaryRow.meanDeltaR_empMinusModel = mean(rEmp - rModel);
            summaryRow.meanEmpiricalR = mean(rEmp);
            summaryRow.meanModelR = mean(rModel);
            summaryRow.meanAbsEmpiricalR = mean(abs(rEmp));
            summaryRow.meanAbsModelR = mean(abs(rModel));
            summaryRow.meanDeltaAbsR_empMinusModel = mean(abs(rEmp) - abs(rModel));
            summaryRow.nPositiveDeltaZ = sum(d > 0);
            summaryRow.nNegativeDeltaZ = sum(d < 0);
            summaryRow.signConsistency = max(summaryRow.nPositiveDeltaZ, summaryRow.nNegativeDeltaZ) / nPatientsUsed;
            summaryRow.meanNvalidParts = mean(nValidParts);
            summaryRow.minNvalidParts = min(nValidParts);
            summaryRow.maxNvalidParts = max(nValidParts);
            summaryRow.errorMessage = "";

            fprintf('  %s vs %s: nPatients = %d, meanDeltaZ = %.4g, t = %.4g, p = %.4g\n', ...
                metabA, metabB, nPatientsUsed, summaryRow.meanDeltaZ_empMinusModel, tValue, pValue);

        else
            summaryRow.status = "failed";
            summaryRow.errorMessage = sprintf('Only %d valid patients; minimum requested is %d.', ...
                nPatientsUsed, opts.minPatientsForGroupTest);
        end

        pairSummaryTable = [pairSummaryTable; summaryRow]; %#ok<AGROW>
    end
end

% Multiple-comparison correction across metabolite pairs.
pairSummaryTable.qValue_FDR = BenjaminiHochberg(pairSummaryTable.pValue);
pairSummaryTable.rejectH0_FDR = pairSummaryTable.qValue_FDR < opts.alpha;

% Sort with most significant pairs first, keeping failed/NaN rows lower.
sortP = pairSummaryTable.pValue;
sortP(~isfinite(sortP)) = Inf;
[~, order] = sort(sortP, 'ascend');
pairSummaryTable = pairSummaryTable(order, :);

pairOutputs = struct();
pairOutputs.method = "Pairwise Fisher-z empirical-vs-LCModel/deGraaf correlation t-test";
pairOutputs.nullHypothesis = "For each metabolite pair, mean over patients of atanh(r_empirical) - atanh(r_model) equals zero.";
pairOutputs.options = opts;
pairOutputs.patientIDs = patientIDs;
pairOutputs.metabolites = metabs;
pairOutputs.crlbMajorityTable = crlbMajorityTable;
pairOutputs.pairSummaryTable = pairSummaryTable;
pairOutputs.patientPairTable = patientPairTable;

if opts.exportResults
    if ~exist(opts.outputDir, 'dir')
        mkdir(opts.outputDir);
    end

    writetable(pairSummaryTable, fullfile(opts.outputDir, 'Pairwise_Empirical_vs_Model_Summary.csv'));
    writetable(patientPairTable, fullfile(opts.outputDir, 'Pairwise_Empirical_vs_Model_PatientLevel.csv'));

    fprintf('\nExported pairwise correlation test results to:\n  %s\n', opts.outputDir);
end

fprintf('\nPairwise empirical-vs-model correlation tests complete.\n');

end

% -------------------------------------------------------------------------
function opts = ParsePairwiseOptions(cfg)

opts = struct();
opts.patientIDs = GetOption(cfg, 'patientIDs', "all");
opts.metabolites = GetOption(cfg, 'metabolites', "all");

opts.excludeSumMetabolites = logical(GetOption(cfg, 'excludeSumMetabolites', true));
opts.sumMetabolites = string(GetOption(cfg, 'sumMetabolites', ...
    ["GPC+PCh", "NAA+NAAG", "Cr+PCr", "Glu+Gln"]));

% CRLB-majority filter:
% keep a metabolite only if, across all selected patients and parts,
% the number of finite CRLB values below the threshold is larger than
% the number of finite CRLB values at/above the threshold.
opts.useCRLBMajorityFilter = logical(GetOption(cfg, 'useCRLBMajorityFilter', true));
opts.crlbMajorityThreshold = double(GetOption(cfg, 'crlbMajorityThreshold', 100));

opts.ignoreZeros = logical(GetOption(cfg, 'ignoreZeros', true));
opts.minValidParts = double(GetOption(cfg, 'minValidParts', 10));
opts.minPatientsForGroupTest = double(GetOption(cfg, 'minPatientsForGroupTest', 3));
opts.alpha = double(GetOption(cfg, 'alpha', 0.05));

opts.exportResults = logical(GetOption(cfg, 'exportResults', false));
opts.outputDir = char(GetOption(cfg, 'outputDir', fullfile(pwd, 'PairwiseEmpiricalVsModelResults')));

end

% -------------------------------------------------------------------------
function value = GetOption(s, fieldName, defaultValue)

if isstruct(s) && isfield(s, fieldName)
    value = s.(fieldName);
else
    value = defaultValue;
end

end

% -------------------------------------------------------------------------
function patientIDs = ResolvePatientIDs(covOutputs, deGraafOutputs, requestedIDs)

if ~isfield(covOutputs, 'patientResultsByID')
    error('covOutputs.patientResultsByID is missing.');
end
if ~isfield(deGraafOutputs, 'patientResultsByID')
    error('deGraafOutputs.patientResultsByID is missing.');
end

covFields = string(fieldnames(covOutputs.patientResultsByID));
modelFields = string(fieldnames(deGraafOutputs.patientResultsByID));
commonFields = intersect(covFields, modelFields, 'stable');

if isempty(commonFields)
    error('No common patient IDs found between covOutputs and deGraafOutputs.');
end

if isnumeric(requestedIDs)
    tmp = strings(numel(requestedIDs), 1);
    for k = 1:numel(requestedIDs)
        tmp(k) = "P" + sprintf('%02d', requestedIDs(k));
    end
    requestedIDs = tmp;
end

requestedIDs = string(requestedIDs);

if isscalar(requestedIDs) && strcmpi(requestedIDs, "all")
    patientIDs = commonFields;
else
    requestedFields = strings(numel(requestedIDs), 1);
    for k = 1:numel(requestedIDs)
        requestedFields(k) = string(matlab.lang.makeValidName(char(requestedIDs(k))));
    end

    keep = ismember(requestedFields, commonFields);
    missing = requestedIDs(~keep);
    if ~isempty(missing)
        warning('Requested patient IDs not found in both outputs: %s', strjoin(cellstr(missing), ', '));
    end

    patientIDs = requestedFields(keep);
end

end

% -------------------------------------------------------------------------
function [metabs, crlbMajorityTable] = ResolveMetabolites(covOutputs, deGraafOutputs, patientIDs, opts)

crlbMajorityTable = table();

if ~(isscalar(string(opts.metabolites)) && strcmpi(string(opts.metabolites), "all"))
    metabs = string(opts.metabolites(:));
else
    % Prefer group-level empirical correlation table, if available.
    if isfield(covOutputs, 'group') && isfield(covOutputs.group, 'meanCorrTable')
        T = covOutputs.group.meanCorrTable;
        if ~isempty(T.Properties.RowNames)
            metabs = string(T.Properties.RowNames);
        else
            metabs = string(T.Properties.VariableNames);
        end
    else
        % Fall back to first patient's partTable variable names.
        firstField = matlab.lang.makeValidName(char(patientIDs(1)));
        T = covOutputs.patientResultsByID.(firstField).partTable;
        metabs = string(T.Properties.VariableNames);
    end

    % Remove non-metabolite bookkeeping columns.
    removeNames = ["part", "patient", "patientID", "division", "file", "filename"];
    keep = true(numel(metabs), 1);
    for k = 1:numel(metabs)
        keep(k) = ~any(strcmpi(metabs(k), removeNames));
    end
    metabs = metabs(keep);
end

if opts.excludeSumMetabolites
    metabs = ApplySumPreferredFiltering(metabs, opts.sumMetabolites);
end

% Apply CRLB-majority filtering before running any pairwise tests.
% This keeps only metabolites with more CRLB values below threshold than
% at/above threshold across all selected patients and parts.
if opts.useCRLBMajorityFilter
    crlbMajorityTable = BuildGlobalCRLBMajorityTable( ...
        deGraafOutputs, patientIDs, opts.crlbMajorityThreshold);

    validCRLBMetabs = crlbMajorityTable.metabolite(crlbMajorityTable.keepByCRLB);
    metabs = KeepMatchingNames(metabs, validCRLBMetabs);

    fprintf('CRLB-majority filter kept %d/%d metabolites using threshold %.1f.\n', ...
        numel(validCRLBMetabs), height(crlbMajorityTable), opts.crlbMajorityThreshold);
end

metabs = metabs(:);

% Keep metabolites that appear in at least one empirical partTable and one
% model table. Per-pair/per-patient validity is still checked later.
keep = false(numel(metabs), 1);
for mIdx = 1:numel(metabs)
    metab = metabs(mIdx);
    for pIdx = 1:numel(patientIDs)
        patientField = matlab.lang.makeValidName(char(patientIDs(pIdx)));
        covPatient = covOutputs.patientResultsByID.(patientField);
        modelPatient = deGraafOutputs.patientResultsByID.(patientField);
        if isfield(covPatient, 'partTable') && HasMatchingTableVariable(covPatient.partTable, metab) && HasMetabInModel(modelPatient, metab)
            keep(mIdx) = true;
            break;
        end
    end
end
metabs = metabs(keep);

if numel(metabs) < 2
    error('Fewer than 2 metabolites available for pairwise testing.');
end

end

% -------------------------------------------------------------------------
function tf = HasMatchingTableVariable(T, label)

varNames = string(T.Properties.VariableNames);
idx = FindMatchingNameIndex(varNames, label);
tf = ~isempty(idx);

end

% -------------------------------------------------------------------------
function tf = HasMetabInModel(modelPatient, metab)

tf = false;

if isfield(modelPatient, 'meanAmplitudeCorrTable')
    T = modelPatient.meanAmplitudeCorrTable;
    labels = GetSquareTableLabels(T);
    vars = string(T.Properties.VariableNames);
    tf = ~isempty(FindMatchingNameIndex(labels, metab)) && ~isempty(FindMatchingNameIndex(vars, metab));
    return;
end

if isfield(modelPatient, 'meanAmplitudeCovTable')
    T = modelPatient.meanAmplitudeCovTable;
    labels = GetSquareTableLabels(T);
    vars = string(T.Properties.VariableNames);
    tf = ~isempty(FindMatchingNameIndex(labels, metab)) && ~isempty(FindMatchingNameIndex(vars, metab));
end

end


% -------------------------------------------------------------------------
function crlbMajorityTable = BuildGlobalCRLBMajorityTable(deGraafOutputs, patientIDs, threshold)

% Build one reliability row per metabolite using all selected patients and
% all available Division_1 parts already stored in deGraafOutputs.
%
% Keep rule:
%   keepByCRLB = nCRLBUnderThreshold > nCRLBOverOrEqualThreshold

if ~isfield(deGraafOutputs, 'patientResultsByID')
    error('deGraafOutputs.patientResultsByID is missing.');
end

allMetabs = strings(0, 1);
allCRLB = cell(0, 1);

for pIdx = 1:numel(patientIDs)

    patientField = matlab.lang.makeValidName(char(patientIDs(pIdx)));

    if ~isfield(deGraafOutputs.patientResultsByID, patientField)
        continue;
    end

    patientResult = deGraafOutputs.patientResultsByID.(patientField);

    if ~isfield(patientResult, 'partCRLB') || ~isfield(patientResult, 'metabList')
        error(['Patient %s does not contain partCRLB/metabList. ', ...
               'Rerun DeGraafAmplitudeCorrelationByPatient with partCRLB saved.'], patientField);
    end

    metabs = string(patientResult.metabList(:));
    partCRLB = double(patientResult.partCRLB);

    for mIdx = 1:numel(metabs)
        existingIdx = FindMatchingNameIndex(allMetabs, metabs(mIdx));

        vals = partCRLB(:, mIdx);
        vals = vals(isfinite(vals));

        if isempty(existingIdx)
            allMetabs(end + 1, 1) = metabs(mIdx); %#ok<AGROW>
            allCRLB{end + 1, 1} = vals(:); %#ok<AGROW>
        else
            allCRLB{existingIdx} = [allCRLB{existingIdx}; vals(:)]; %#ok<AGROW>
        end
    end
end

crlbMajorityTable = table();

for mIdx = 1:numel(allMetabs)

    vals = allCRLB{mIdx};
    vals = vals(isfinite(vals));

    nFinite = numel(vals);
    nUnder = sum(vals < threshold);
    nOverOrEqual = sum(vals >= threshold);

    row = table();
    row.metabolite = allMetabs(mIdx);
    row.nFiniteCRLB = nFinite;
    row.nCRLBUnderThreshold = nUnder;
    row.nCRLBOverOrEqualThreshold = nOverOrEqual;
    row.percentUnderThreshold = 100 * nUnder / max(nFinite, 1);
    row.threshold = threshold;
    row.keepByCRLB = nUnder > nOverOrEqual;

    crlbMajorityTable = [crlbMajorityTable; row]; %#ok<AGROW>
end

end

% -------------------------------------------------------------------------
function metabsOut = ApplySumPreferredFiltering(metabsIn, sumMetabolites)

% Sum-preferred rule:
% If a summed metabolite exists, keep the sum and remove its components.

metabsOut = string(metabsIn(:));
sumMetabolites = string(sumMetabolites(:));

for sIdx = 1:numel(sumMetabolites)

    sumName = string(sumMetabolites(sIdx));

    % Only remove components if the summed metabolite actually exists.
    sumExists = ~isempty(FindMatchingNameIndex(metabsOut, sumName));

    if ~sumExists
        continue;
    end

    components = string(split(sumName, '+'));
    remove = false(numel(metabsOut), 1);

    for mIdx = 1:numel(metabsOut)
        for cIdx = 1:numel(components)
            if NamesMatch(metabsOut(mIdx), components(cIdx))
                remove(mIdx) = true;
            end
        end
    end

    metabsOut = metabsOut(~remove);
end

end

% -------------------------------------------------------------------------
function metabsOut = KeepMatchingNames(metabsIn, allowedMetabs)

metabsIn = string(metabsIn(:));
allowedMetabs = string(allowedMetabs(:));

keep = false(numel(metabsIn), 1);

for mIdx = 1:numel(metabsIn)
    keep(mIdx) = ~isempty(FindMatchingNameIndex(allowedMetabs, metabsIn(mIdx)));
end

metabsOut = metabsIn(keep);

end

% -------------------------------------------------------------------------
function tf = NamesMatch(a, b)

a = string(a);
b = string(b);

if strcmp(a, b)
    tf = true;
    return;
end

aValid = string(matlab.lang.makeValidName(char(a)));
bValid = string(matlab.lang.makeValidName(char(b)));

tf = strcmp(aValid, bValid);

end

% -------------------------------------------------------------------------
function row = MakeEmptyPatientPairRow(patientID, metabA, metabB)

row = table( ...
    string(patientID), ...
    string(metabA), ...
    string(metabB), ...
    "not_run", ...
    NaN, ...
    NaN, ...
    NaN, ...
    NaN, ...
    NaN, ...
    NaN, ...
    NaN, ...
    "", ...
    'VariableNames', { ...
        'patientID', ...
        'metaboliteA', ...
        'metaboliteB', ...
        'status', ...
        'nValidParts', ...
        'rEmpirical', ...
        'rModel', ...
        'zEmpirical', ...
        'zModel', ...
        'deltaZ_empMinusModel', ...
        'deltaR_empMinusModel', ...
        'errorMessage'});

end

% -------------------------------------------------------------------------
function row = MakeEmptyPairSummaryRow(metabA, metabB)

% Build the empty row by named assignments rather than by passing a long
% value list to table(...). This avoids MATLAB errors if a placeholder value
% is accidentally missed.
row = table();

row.metaboliteA = string(metabA);
row.metaboliteB = string(metabB);
row.status = "not_run";

row.nPatientsUsed = NaN;
row.df = NaN;
row.meanDeltaZ_empMinusModel = NaN;
row.sdDeltaZ_empMinusModel = NaN;
row.tValue = NaN;
row.pValue = NaN;
row.rejectH0 = false;

row.groupEmpiricalR = NaN;
row.groupModelR = NaN;
row.meanDeltaR_empMinusModel = NaN;
row.meanEmpiricalR = NaN;
row.meanModelR = NaN;
row.meanAbsEmpiricalR = NaN;
row.meanAbsModelR = NaN;
row.meanDeltaAbsR_empMinusModel = NaN;

row.nPositiveDeltaZ = NaN;
row.nNegativeDeltaZ = NaN;
row.signConsistency = NaN;

row.meanNvalidParts = NaN;
row.minNvalidParts = NaN;
row.maxNvalidParts = NaN;

row.errorMessage = "";

end

% -------------------------------------------------------------------------
function values = GetTableColumnByMetab(T, metab)

varNames = string(T.Properties.VariableNames);
idx = FindMatchingNameIndex(varNames, metab);

if isempty(idx)
    error('Could not find metabolite %s in table.', metab);
end

values = double(T{:, idx});
values = values(:);

end

% -------------------------------------------------------------------------
function rModel = GetModelCorrelation(modelPatient, metabA, metabB)
    
    if isfield(modelPatient, 'meanAmplitudeCorrTable')
        T = modelPatient.meanAmplitudeCorrTable;
        rModel = GetSquareTableValue(T, metabA, metabB);
        return;
    end
    
    if isfield(modelPatient, 'meanAmplitudeCovTable')
        T = modelPatient.meanAmplitudeCovTable;
        labels = GetSquareTableLabels(T);
        M = ExtractSquareMatrixFromTable(T, labels);
        R = CovToCorr(M);
    
        idxA = FindMatchingNameIndex(labels, metabA);
        idxB = FindMatchingNameIndex(labels, metabB);
    
        if isempty(idxA) || isempty(idxB)
            error('Could not find model covariance labels for %s and %s.', metabA, metabB);
        end
    
        rModel = R(idxA, idxB);
        return;
    end
    
    error('Model patient result has neither meanAmplitudeCorrTable nor meanAmplitudeCovTable.');

end

% -------------------------------------------------------------------------
function value = GetSquareTableValue(T, rowLabel, colLabel)

rowLabels = GetSquareTableLabels(T);
colLabels = string(T.Properties.VariableNames);

rowIdx = FindMatchingNameIndex(rowLabels, rowLabel);
colIdx = FindMatchingNameIndex(colLabels, colLabel);

if isempty(rowIdx)
    error('Could not find row %s in square table.', rowLabel);
end
if isempty(colIdx)
    error('Could not find column %s in square table.', colLabel);
end

value = double(T{rowIdx, colIdx});

end

% -------------------------------------------------------------------------
function labels = GetSquareTableLabels(T)

if ~istable(T)
    error('Expected a MATLAB table.');
end

if ~isempty(T.Properties.RowNames)
    labels = string(T.Properties.RowNames);
else
    labels = string(T.Properties.VariableNames);
end

labels = labels(:);

end

% -------------------------------------------------------------------------
function M = ExtractSquareMatrixFromTable(T, labels)

labels = string(labels(:));
rowLabels = GetSquareTableLabels(T);
varNames = string(T.Properties.VariableNames);

rowIdx = nan(numel(labels), 1);
colIdx = nan(numel(labels), 1);

for k = 1:numel(labels)
    r = FindMatchingNameIndex(rowLabels, labels(k));
    c = FindMatchingNameIndex(varNames, labels(k));

    if isempty(r)
        error('Could not find row %s in square table.', labels(k));
    end
    if isempty(c)
        error('Could not find column %s in square table.', labels(k));
    end

    rowIdx(k) = r;
    colIdx(k) = c;
end

M = double(T{rowIdx, colIdx});
M = (M + M.') / 2;

end

% -------------------------------------------------------------------------
function idx = FindMatchingNameIndex(names, label)

names = string(names(:));
label = string(label);

idx = find(strcmp(names, label), 1, 'first');
if ~isempty(idx)
    return;
end

labelValid = string(matlab.lang.makeValidName(char(label)));
idx = find(strcmp(names, labelValid), 1, 'first');
if ~isempty(idx)
    return;
end

namesValid = strings(numel(names), 1);
for k = 1:numel(names)
    namesValid(k) = string(matlab.lang.makeValidName(char(names(k))));
end
idx = find(strcmp(namesValid, labelValid), 1, 'first');

end

% -------------------------------------------------------------------------
function r = PearsonR(x, y)

x = double(x(:));
y = double(y(:));

valid = isfinite(x) & isfinite(y);
x = x(valid);
y = y(valid);

if numel(x) < 3
    r = NaN;
    return;
end

x = x - mean(x);
y = y - mean(y);

denom = sqrt(sum(x.^2) * sum(y.^2));

if denom <= 0
    r = NaN;
else
    r = sum(x .* y) / denom;
end

end

% -------------------------------------------------------------------------
function r = ClampCorrelation(r)

r = max(min(r, 0.999999), -0.999999);

end

% -------------------------------------------------------------------------
function [tValue, pValue, df] = OneSampleTTestAgainstZero(x)

x = double(x(:));
x = x(isfinite(x));

n = numel(x);
df = n - 1;

if n < 3
    tValue = NaN;
    pValue = NaN;
    return;
end

sx = std(x, 0);
if sx <= 0 || ~isfinite(sx)
    tValue = NaN;
    pValue = NaN;
    return;
end

mx = mean(x);
tValue = mx / (sx / sqrt(n));

try
    pValue = 2 * tcdf(-abs(tValue), df);
catch
    pValue = NaN;
    warning('tcdf is unavailable. pValue was set to NaN.');
end

end

% -------------------------------------------------------------------------
function q = BenjaminiHochberg(p)

p = p(:);
q = nan(size(p));

valid = isfinite(p);
pValid = p(valid);

if isempty(pValid)
    return;
end

[pSorted, sortIdx] = sort(pValid, 'ascend');
m = numel(pSorted);

qSorted = pSorted .* m ./ (1:m)';

for k = m-1:-1:1
    qSorted(k) = min(qSorted(k), qSorted(k+1));
end

qSorted = min(qSorted, 1);

qValid = nan(size(pValid));
qValid(sortIdx) = qSorted;

q(valid) = qValid;

end

% -------------------------------------------------------------------------
function R = CovToCorr(Sigma)

Sigma = double(Sigma);
sd = sqrt(diag(Sigma));
R = Sigma ./ (sd * sd.');
R(1:size(R,1)+1:end) = 1;

end
