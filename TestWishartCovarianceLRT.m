function lrtOutputs = TestWishartCovarianceLRT(inputA, varargin)
% TestWishartCovarianceLRT
%
% Global covariance-matrix likelihood-ratio test for repeated MRS fits.
%
% Null hypothesis for each patient:
%   H0: Sigma_true = Sigma_model
%
% Test statistic:
%   T = nu * [ trace(Sigma_model^-1 * S_emp)
%              - log(det(Sigma_model^-1 * S_emp))
%              - p ]
%
% Approximate null distribution:
%   T ~ chi-square(df)
%   df = p*(p + 1)/2
%
% Metabolite-selection modes:
%   lrtCfg.metaboliteSelectionMode = "perPatientLargestValid";
%       Mode A. Each patient uses the largest valid metabolite subset available
%       for that patient. Patients may therefore use different panels.
%
%   lrtCfg.metaboliteSelectionMode = "fixed";
%       Mode B. Use a predefined metabolite list from lrtCfg.metabolites.
%       All patients are tested with the same requested panel; patients that
%       cannot support it fail. Alias: "predefined".
%
%   lrtCfg.metaboliteSelectionMode = "largestCommon";
%       Mode C. First find the largest valid panel per patient, then take the
%       common intersection across eligible patients, then rerun the LRT using
%       that one common panel for every included patient.

[analysisData, lrtCfg] = ResolveCentralizedInputs(inputA, varargin{:});
[covOutputs, deGraafOutputs] = ExtractCanonicalOutputs(analysisData);
opts = ParseLRTOptions(lrtCfg);

viewName = char(string(GetOption(lrtCfg, 'filterView', "default")));
if ~isfield(analysisData.wishart.views, viewName)
    error('Unknown centralized Wishart filter view: %s.', viewName);
end
filterView = analysisData.wishart.views.(viewName);
patientIDs = filterView.patientIDs;

% Scientific filtering was completed by ApplyAnalysisFilters. These values
% prevent the legacy selection helpers below from applying it a second time.
opts.patientIDs = patientIDs;
opts.metabolites = filterView.metabolites;
opts.minValidParts = filterView.minValidParts;
opts.crlbMajorityTable = analysisData.crlbMajorityTable;
opts.crlbValidMetabolites = analysisData.crlbMajorityTable.metabolite( ...
    analysisData.crlbMajorityTable.keepByCRLB);
crlbMajorityTable = analysisData.crlbMajorityTable;

panelDiscoveryTable = table();
commonPanel = strings(0, 1);
patientIDsToRun = patientIDs;
optsRun = opts;

if strcmpi(opts.metaboliteSelectionMode, "largestCommon")
    fprintf('\nDiscovering largest common valid metabolite panel...\n');

    [commonPanel, panelDiscoveryTable, eligiblePatientIDs] = DiscoverLargestCommonPanel( ...
        covOutputs, deGraafOutputs, patientIDs, opts);

    optsRun = opts;
    optsRun.metaboliteSelectionMode = "fixed";
    optsRun.metabolites = commonPanel;

    if opts.runOnlyCommonPanelPatients
        patientIDsToRun = eligiblePatientIDs;
    else
        patientIDsToRun = patientIDs;
    end

    fprintf('\nCommon metabolite panel used for final LRT (%d metabolites):\n', numel(commonPanel));
    fprintf('  %s\n', strjoin(cellstr(commonPanel), ', '));
    fprintf('Final LRT will run for %d patient(s).\n', numel(patientIDsToRun));
end

[patientResultsByID, summaryRows, groupTable] = RunPatientsLRT( ...
    covOutputs, deGraafOutputs, patientIDsToRun, optsRun);

lrtOutputs = struct();
lrtOutputs.method = "Wishart covariance likelihood-ratio test";
lrtOutputs.nullHypothesis = "For each patient, Sigma_true equals the LCModel/deGraaf model covariance matrix.";
lrtOutputs.testStatistic = "T = nu * [trace(Sigma0^-1*S) - logdet(Sigma0^-1*S) - p]";
lrtOutputs.dfDefinition = "df = p*(p+1)/2 for full covariance matrix";
lrtOutputs.options = opts;
lrtOutputs.patientIDs = patientIDsToRun;
lrtOutputs.allCandidatePatientIDs = patientIDs;
lrtOutputs.metaboliteSelectionMode = opts.metaboliteSelectionMode;
lrtOutputs.commonPanel = commonPanel;
lrtOutputs.panelDiscoveryTable = panelDiscoveryTable;
lrtOutputs.crlbMajorityTable = crlbMajorityTable;
lrtOutputs.patientResultsByID = patientResultsByID;
lrtOutputs.patientSummaryTable = summaryRows;
lrtOutputs.groupTable = groupTable;

if opts.exportResults
    if ~exist(opts.outputDir, 'dir')
        mkdir(opts.outputDir);
    end

    writetable(summaryRows, fullfile(opts.outputDir, 'Wishart_LRT_PerPatient.csv'));
    writetable(groupTable, fullfile(opts.outputDir, 'Wishart_LRT_Group.csv'));

    if ~isempty(panelDiscoveryTable)
        writetable(panelDiscoveryTable, fullfile(opts.outputDir, 'Wishart_LRT_CommonPanelDiscovery.csv'));
    end

    if ~isempty(commonPanel)
        commonPanelTable = table(commonPanel(:), 'VariableNames', {'metabolite'});
        writetable(commonPanelTable, fullfile(opts.outputDir, 'Wishart_LRT_CommonPanel.csv'));
    end

    fprintf('\nExported Wishart LRT results to:\n  %s\n', opts.outputDir);
end

fprintf('\nWishart covariance LRT complete.\n');
if ~isempty(groupTable) && isfinite(groupTable.T_group(1))
    fprintf('Group result: T = %.4g, df = %d, p = %.4g\n', ...
        groupTable.T_group(1), groupTable.df_group(1), groupTable.pValue_group(1));
end

end

% -------------------------------------------------------------------------
function [patientResultsByID, summaryRows, groupTable] = RunPatientsLRT(covOutputs, deGraafOutputs, patientIDs, opts)

patientResultsByID = struct();
summaryRows = table();

fprintf('\nRunning Wishart covariance LRT for %d patient(s)...\n', numel(patientIDs));

for pIdx = 1:numel(patientIDs)

    patientID = string(patientIDs(pIdx));
    patientField = matlab.lang.makeValidName(char(patientID));

    row = MakeEmptySummaryRow(patientID);

    try
        [partTable, modelCovTable] = GetPatientTables(covOutputs, deGraafOutputs, patientID);

        optsPatient = opts;
        if strcmpi(opts.metaboliteSelectionMode, "perPatientLargestValid")
            candidateMetabs = SelectCandidateMetabolites(partTable, modelCovTable, opts);
            validPanel = FindLargestValidPanelForPatient(partTable, modelCovTable, candidateMetabs, opts);
            optsPatient.metaboliteSelectionMode = "fixed";
            optsPatient.metabolites = validPanel;
        end

        patientResult = RunOnePatientLRT(partTable, modelCovTable, patientID, optsPatient);

        patientResultsByID.(patientField) = patientResult;
        row = patientResult.summaryRow;

        fprintf('  %s: T = %.4g, df = %d, p = %.4g, pMetabs = %d, Nvalid = %d\n', ...
            patientID, row.T, row.df, row.pValue, row.nMetabolites, row.nValidParts);

    catch ME
        warning('Wishart covariance LRT failed for patient %s: %s', patientID, ME.message);
        row.status = "failed";
        row.errorMessage = string(ME.message);
        patientResultsByID.(patientField).summaryRow = row;
    end

    summaryRows = [summaryRows; row]; %#ok<AGROW>
end

validForGroup = strcmp(summaryRows.status, "ok") & isfinite(summaryRows.T) & isfinite(summaryRows.df);

if any(validForGroup)
    groupT = sum(summaryRows.T(validForGroup));
    groupDf = sum(summaryRows.df(validForGroup));
    groupP = ChiSquareUpperTail(groupT, groupDf);
    groupReject = groupP < opts.alpha;
else
    groupT = NaN;
    groupDf = NaN;
    groupP = NaN;
    groupReject = false;
end

groupTable = table( ...
    sum(validForGroup), ...
    groupT, ...
    groupDf, ...
    groupP, ...
    groupReject, ...
    opts.alpha, ...
    'VariableNames', {'nPatientsUsed', 'T_group', 'df_group', 'pValue_group', 'rejectH0_group', 'alpha'});

end

% -------------------------------------------------------------------------
function [commonPanel, panelDiscoveryTable, eligiblePatientIDs] = DiscoverLargestCommonPanel(covOutputs, deGraafOutputs, patientIDs, opts)

panelDiscoveryTable = table();
eligiblePatientIDs = strings(0, 1);
commonPanel = strings(0, 1);
firstEligible = true;

for pIdx = 1:numel(patientIDs)

    patientID = string(patientIDs(pIdx));
    row = MakeEmptyDiscoveryRow(patientID);

    try
        [partTable, modelCovTable] = GetPatientTables(covOutputs, deGraafOutputs, patientID);

        optsDiscovery = opts;
        optsDiscovery.metabolites = opts.metabolites;

        candidateMetabs = SelectCandidateMetabolites(partTable, modelCovTable, optsDiscovery);
        row.nCandidateMetabolites = numel(candidateMetabs);

        validPanel = FindLargestValidPanelForPatient(partTable, modelCovTable, candidateMetabs, opts);

        row.status = "ok";
        row.isEligible = true;
        row.nValidMetabolites = numel(validPanel);
        row.validMetabolites = string(strjoin(cellstr(validPanel(:)), ', '));
        row.errorMessage = "";

        eligiblePatientIDs(end + 1, 1) = string(matlab.lang.makeValidName(char(patientID))); %#ok<AGROW>

        if firstEligible
            commonPanel = validPanel(:);
            firstEligible = false;
        else
            keep = ismember(commonPanel, validPanel);
            commonPanel = commonPanel(keep);
        end

        fprintf('  %s eligible panel (%d): %s\n', ...
            patientID, numel(validPanel), strjoin(cellstr(validPanel), ', '));

    catch ME
        row.status = "failed";
        row.isEligible = false;
        row.errorMessage = string(ME.message);
        fprintf('  %s not eligible for panel discovery: %s\n', patientID, ME.message);
    end

    panelDiscoveryTable = [panelDiscoveryTable; row]; %#ok<AGROW>
end

if isempty(eligiblePatientIDs)
    error('No eligible patients found during common-panel discovery.');
end

commonPanel = commonPanel(:);

if numel(commonPanel) < opts.minMetabolites
    error('Common panel contains only %d metabolite(s); minimum requested is %d.', ...
        numel(commonPanel), opts.minMetabolites);
end

end

% -------------------------------------------------------------------------
function validPanel = FindLargestValidPanelForPatient(partTable, modelCovTable, candidateMetabs, opts)

candidateMetabs = string(candidateMetabs(:));

if isempty(candidateMetabs)
    error('No candidate metabolites available.');
end

currentMetabs = candidateMetabs;
lastError = "";

while numel(currentMetabs) >= opts.minMetabolites

    try
        [X, dataLabels] = ExtractDataMatrixFromPartTable(partTable, currentMetabs);
        Sigma0 = ExtractSquareMatrixFromTable(modelCovTable, dataLabels);

        % Remove metabolites whose model covariance row/column is not usable.
        goodModel = all(isfinite(Sigma0), 1)' & all(isfinite(Sigma0), 2) & ...
                    isfinite(diag(Sigma0)) & diag(Sigma0) > 0;
        X = X(:, goodModel);
        Sigma0 = Sigma0(goodModel, goodModel);
        dataLabels = dataLabels(goodModel);

        if numel(dataLabels) < opts.minMetabolites
            error('Fewer than %d metabolites remained after model-covariance filtering.', opts.minMetabolites);
        end

        validRows = all(isfinite(X), 2);
        Xvalid = X(validRows, :);
        Nvalid = size(Xvalid, 1);

        if Nvalid < opts.minValidParts
            % Try removing the metabolite with the most missing values.
            missingCounts = sum(~isfinite(X), 1);
            [maxMissing, worstIdx] = max(missingCounts);
            if maxMissing <= 0
                error('Only %d valid repeated parts remained; removing metabolites cannot improve this.', Nvalid);
            end
            currentMetabs(strcmp(currentMetabs, dataLabels(worstIdx))) = [];
            lastError = sprintf('Nvalid = %d < minValidParts = %d.', Nvalid, opts.minValidParts);
            continue;
        end

        empVar = var(Xvalid, 0, 1);
        goodEmpirical = isfinite(empVar(:)) & empVar(:) > 0;
        Xvalid = Xvalid(:, goodEmpirical);
        Sigma0 = Sigma0(goodEmpirical, goodEmpirical);
        dataLabels = dataLabels(goodEmpirical);

        p = numel(dataLabels);

        if p < opts.minMetabolites
            error('Fewer than %d metabolites remained after empirical variance filtering.', opts.minMetabolites);
        end

        if Nvalid <= p
            % Remove one metabolite to make covariance full rank.
            missingCounts = sum(~isfinite(X), 1);
            missingCounts = missingCounts(goodEmpirical);
            if isempty(missingCounts)
                worstIdx = p;
            else
                [~, worstIdx] = max(missingCounts);
                if all(missingCounts == 0)
                    worstIdx = p;
                end
            end
            currentMetabs(strcmp(currentMetabs, dataLabels(worstIdx))) = [];
            lastError = sprintf('Need Nvalid > p. Got Nvalid = %d and p = %d.', Nvalid, p);
            continue;
        end

        Semp = ForceSymmetric(cov(Xvalid, 0));
        Sigma0 = ForceSymmetric(Sigma0);

        % Check positive definiteness / ridge compatibility.
        MakePositiveDefiniteIfNeeded(Semp, opts, 'empirical covariance');
        MakePositiveDefiniteIfNeeded(Sigma0, opts, 'model covariance');

        validPanel = dataLabels(:);
        return;

    catch ME
        lastError = string(ME.message);

        % If current panel cannot even be evaluated, remove one metabolite and retry.
        if numel(currentMetabs) <= opts.minMetabolites
            break;
        end
        currentMetabs(end) = [];
    end
end

error('Could not find a valid panel for this patient. Last reason: %s', lastError);

end

% -------------------------------------------------------------------------
function [partTable, modelCovTable] = GetPatientTables(covOutputs, deGraafOutputs, patientID)

patientField = matlab.lang.makeValidName(char(patientID));

if ~isfield(covOutputs, 'patientResultsByID') || ~isfield(covOutputs.patientResultsByID, patientField)
    error('Missing covOutputs.patientResultsByID.%s.', patientField);
end

if ~isfield(deGraafOutputs, 'patientResultsByID') || ~isfield(deGraafOutputs.patientResultsByID, patientField)
    error('Missing deGraafOutputs.patientResultsByID.%s.', patientField);
end

covPatient = covOutputs.patientResultsByID.(patientField);
modelPatient = deGraafOutputs.patientResultsByID.(patientField);

if ~isfield(covPatient, 'partTable')
    error('Missing partTable for patient %s.', patientID);
end

if ~isfield(modelPatient, 'meanAmplitudeCovTable')
    error('Missing meanAmplitudeCovTable for patient %s.', patientID);
end

partTable = covPatient.partTable;
modelCovTable = modelPatient.meanAmplitudeCovTable;

end

% -------------------------------------------------------------------------
function opts = ParseLRTOptions(cfg)

opts = struct();
modeDefault = "perPatientLargestValid";
if isstruct(cfg) && isfield(cfg, 'useLargestCommonSubset') && logical(cfg.useLargestCommonSubset)
    modeDefault = "largestCommon";
end
opts.metaboliteSelectionMode = string(GetOption(cfg, 'metaboliteSelectionMode', modeDefault));

validModes = ["perPatientLargestValid", "perPatient", "a", ...
              "fixed", "predefined", "b", ...
              "largestCommon", "common", "c"];
if ~any(strcmpi(opts.metaboliteSelectionMode, validModes))
    error(['Unknown metaboliteSelectionMode: %s. Use "perPatientLargestValid" ', ...
           '(or "perPatient"/"a"), "fixed" (or "predefined"/"b"), ', ...
           'or "largestCommon" (or "common"/"c").'], opts.metaboliteSelectionMode);
end

% Normalize spelling/case.
if strcmpi(opts.metaboliteSelectionMode, "perPatientLargestValid") || ...
   strcmpi(opts.metaboliteSelectionMode, "perPatient") || ...
   strcmpi(opts.metaboliteSelectionMode, "a")
    opts.metaboliteSelectionMode = "perPatientLargestValid";
elseif strcmpi(opts.metaboliteSelectionMode, "fixed") || ...
       strcmpi(opts.metaboliteSelectionMode, "predefined") || ...
       strcmpi(opts.metaboliteSelectionMode, "b")
    opts.metaboliteSelectionMode = "fixed";
elseif strcmpi(opts.metaboliteSelectionMode, "largestCommon") || ...
       strcmpi(opts.metaboliteSelectionMode, "common") || ...
       strcmpi(opts.metaboliteSelectionMode, "c")
    opts.metaboliteSelectionMode = "largestCommon";
end

opts.runOnlyCommonPanelPatients = logical(GetOption(cfg, 'runOnlyCommonPanelPatients', true));
opts.minMetabolites = double(GetOption(cfg, 'minMetabolites', 2));

opts.crlbMajorityTable = table();
opts.crlbValidMetabolites = strings(0, 1);

opts.alpha = double(GetOption(cfg, 'alpha', 0.05));

% Numerical ridge is only meant to fix tiny numerical non-positive-definite
% problems. If a matrix needs a large ridge, inspect the metabolites used.
opts.applyNumericalRidge = logical(GetOption(cfg, 'applyNumericalRidge', true));
opts.ridgeScale = double(GetOption(cfg, 'ridgeScale', 1e-8));
opts.maxRidgeSteps = double(GetOption(cfg, 'maxRidgeSteps', 10));

opts.exportResults = logical(GetOption(cfg, 'exportResults', false));
opts.outputDir = char(GetOption(cfg, 'outputDir', fullfile(pwd, 'WishartLRTResults')));

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
function patientResult = RunOnePatientLRT(partTable, modelCovTable, patientID, opts)

candidateMetabs = SelectCandidateMetabolites(partTable, modelCovTable, opts);

if isempty(candidateMetabs)
    error('No candidate metabolites remained after matching partTable and model covariance table.');
end

[X, dataLabels] = ExtractDataMatrixFromPartTable(partTable, candidateMetabs);
Sigma0 = ExtractSquareMatrixFromTable(modelCovTable, dataLabels);

% Remove metabolites with non-finite model entries.
goodModel = all(isfinite(Sigma0), 1)' & all(isfinite(Sigma0), 2) & isfinite(diag(Sigma0)) & diag(Sigma0) > 0;
X = X(:, goodModel);
Sigma0 = Sigma0(goodModel, goodModel);
dataLabels = dataLabels(goodModel);

if isempty(dataLabels)
    error('No metabolites remained after removing non-finite model covariance entries.');
end

% Listwise-complete rows are required for a single valid covariance matrix.
validRows = all(isfinite(X), 2);
Xvalid = X(validRows, :);
Nvalid = size(Xvalid, 1);

if Nvalid < opts.minValidParts
    error('Only %d valid repeated parts remained; minimum requested is %d.', Nvalid, opts.minValidParts);
end

% Remove zero-variance empirical columns.
empVar = var(Xvalid, 0, 1);
goodEmpirical = isfinite(empVar(:)) & empVar(:) > 0;
Xvalid = Xvalid(:, goodEmpirical);
Sigma0 = Sigma0(goodEmpirical, goodEmpirical);
dataLabels = dataLabels(goodEmpirical);

p = numel(dataLabels);

if p < opts.minMetabolites
    error('Fewer than %d metabolites remained; matrix test cannot be run.', opts.minMetabolites);
end

if Nvalid <= p
    error(['Need Nvalid > p for a full-rank empirical covariance matrix. ' ...
           'Got Nvalid = %d and p = %d. Reduce metabolites or handle missing values.'], Nvalid, p);
end

Semp = cov(Xvalid, 0);  % unbiased covariance: denominator Nvalid - 1
Semp = ForceSymmetric(Semp);
Sigma0 = ForceSymmetric(Sigma0);

[SempPD, ridgeS, cholS] = MakePositiveDefiniteIfNeeded(Semp, opts, 'empirical covariance');
[Sigma0PD, ridgeSigma0, cholSigma0] = MakePositiveDefiniteIfNeeded(Sigma0, opts, 'model covariance');

nu = Nvalid - 1;

traceTerm = trace(Sigma0PD \ SempPD);
logDetS = 2 * sum(log(diag(cholS)));
logDetSigma0 = 2 * sum(log(diag(cholSigma0)));
logDetRatio = logDetS - logDetSigma0;

T = nu * (traceTerm - logDetRatio - p);

% T should be non-negative. Very small negative values can occur numerically.
if T < 0 && abs(T) < 1e-8
    T = 0;
end

df = p * (p + 1) / 2;
pValue = ChiSquareUpperTail(T, df);
rejectH0 = pValue < opts.alpha;

summaryRow = MakeEmptySummaryRow(patientID);
summaryRow.status = "ok";
summaryRow.nValidParts = Nvalid;
summaryRow.nu = nu;
summaryRow.nMetabolites = p;
summaryRow.df = df;
summaryRow.T = T;
summaryRow.pValue = pValue;
summaryRow.rejectH0 = rejectH0;
summaryRow.traceTerm = traceTerm;
summaryRow.logDetRatio = logDetRatio;
summaryRow.ridgeEmpirical = ridgeS;
summaryRow.ridgeModel = ridgeSigma0;
summaryRow.metabolitesUsed = string(strjoin(cellstr(dataLabels(:)), ', '));
summaryRow.nMetabolitesRemoved = numel(candidateMetabs) - p;
summaryRow.errorMessage = "";

patientResult = struct();
patientResult.patientID = patientID;
patientResult.status = "ok";
patientResult.summaryRow = summaryRow;
patientResult.metabolitesUsed = dataLabels;
patientResult.validRows = validRows;
patientResult.Xvalid = Xvalid;
patientResult.S_empirical = Semp;
patientResult.S_empirical_PD_usedForTest = SempPD;
patientResult.Sigma_model = Sigma0;
patientResult.Sigma_model_PD_usedForTest = Sigma0PD;
patientResult.T = T;
patientResult.df = df;
patientResult.pValue = pValue;
patientResult.rejectH0 = rejectH0;
patientResult.nu = nu;
patientResult.nValidParts = Nvalid;
patientResult.traceTerm = traceTerm;
patientResult.logDetRatio = logDetRatio;
patientResult.ridgeEmpirical = ridgeS;
patientResult.ridgeModel = ridgeSigma0;

end

% -------------------------------------------------------------------------
function candidateMetabs = SelectCandidateMetabolites(partTable, modelCovTable, opts)

modelLabels = GetSquareTableLabels(modelCovTable);
if isempty(modelLabels)
    error('Could not determine labels from model covariance table.');
end

candidateMetabs = strings(0, 1);
for k = 1:numel(modelLabels)
    lab = string(modelLabels(k));
    if HasMatchingTableVariable(partTable, lab) && HasMatchingSquareTableColumn(modelCovTable, lab)
        candidateMetabs(end + 1, 1) = lab; %#ok<AGROW>
    end
end

if ~(isscalar(string(opts.metabolites)) && strcmpi(string(opts.metabolites), "all"))
    requested = string(opts.metabolites(:));
    keep = false(numel(candidateMetabs), 1);
    requestedValid = strings(numel(requested), 1);
    for k = 1:numel(requested)
        requestedValid(k) = string(matlab.lang.makeValidName(char(requested(k))));
    end

    for k = 1:numel(candidateMetabs)
        candidateValid = string(matlab.lang.makeValidName(char(candidateMetabs(k))));
        keep(k) = any(strcmp(candidateMetabs(k), requested)) || any(strcmp(candidateValid, requestedValid));
    end
    candidateMetabs = candidateMetabs(keep);
end

candidateMetabs = candidateMetabs(:);

end


% -------------------------------------------------------------------------
function row = MakeEmptySummaryRow(patientID)

row = table( ...
    string(patientID), ...
    "not_run", ...
    NaN, ...
    NaN, ...
    NaN, ...
    NaN, ...
    NaN, ...
    NaN, ...
    false, ...
    NaN, ...
    NaN, ...
    NaN, ...
    NaN, ...
    NaN, ...
    "", ...
    "", ...
    'VariableNames', { ...
        'patientID', ...
        'status', ...
        'nValidParts', ...
        'nu', ...
        'nMetabolites', ...
        'df', ...
        'T', ...
        'pValue', ...
        'rejectH0', ...
        'traceTerm', ...
        'logDetRatio', ...
        'ridgeEmpirical', ...
        'ridgeModel', ...
        'nMetabolitesRemoved', ...
        'metabolitesUsed', ...
        'errorMessage'});

end

% -------------------------------------------------------------------------
function row = MakeEmptyDiscoveryRow(patientID)

row = table( ...
    string(patientID), ...
    "not_run", ...
    false, ...
    NaN, ...
    NaN, ...
    "", ...
    "", ...
    'VariableNames', { ...
        'patientID', ...
        'status', ...
        'isEligible', ...
        'nCandidateMetabolites', ...
        'nValidMetabolites', ...
        'validMetabolites', ...
        'errorMessage'});

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
function tf = HasMatchingTableVariable(T, label)

varNames = string(T.Properties.VariableNames);
idx = FindMatchingNameIndex(varNames, label);
tf = ~isempty(idx);

end

% -------------------------------------------------------------------------
function tf = HasMatchingSquareTableColumn(T, label)

varNames = string(T.Properties.VariableNames);
idx = FindMatchingNameIndex(varNames, label);
tf = ~isempty(idx);

end

% -------------------------------------------------------------------------
function [X, labels] = ExtractDataMatrixFromPartTable(partTable, labels)

labels = string(labels(:));
X = nan(height(partTable), numel(labels));

for k = 1:numel(labels)
    varNames = string(partTable.Properties.VariableNames);
    idx = FindMatchingNameIndex(varNames, labels(k));
    if isempty(idx)
        error('Could not find metabolite %s in partTable.', labels(k));
    end

    values = partTable{:, idx};
    values = double(values);

    X(:, k) = values;
end

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

M = T{rowIdx, colIdx};
M = double(M);
M = ForceSymmetric(M);

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
function A = ForceSymmetric(A)

A = (A + A.') / 2;

end

% -------------------------------------------------------------------------
function [Aout, ridgeUsed, cholA] = MakePositiveDefiniteIfNeeded(A, opts, matrixName)

A = ForceSymmetric(A);
ridgeUsed = 0;

if any(~isfinite(A(:)))
    error('%s contains non-finite values.', matrixName);
end

[cholA, flag] = chol(A, 'lower');
if flag == 0
    Aout = A;
    return;
end

if ~opts.applyNumericalRidge
    error('%s is not positive definite. Try fewer metabolites or exclude sums.', matrixName);
end

diagVals = abs(diag(A));
diagVals = diagVals(isfinite(diagVals) & diagVals > 0);
if isempty(diagVals)
    baseScale = 1;
else
    baseScale = median(diagVals);
end

for step = 0:opts.maxRidgeSteps
    ridge = opts.ridgeScale * (10^step) * max(baseScale, eps);
    Atry = A + ridge * eye(size(A));
    [cholTry, flag] = chol(Atry, 'lower');
    if flag == 0
        Aout = Atry;
        ridgeUsed = ridge;
        cholA = cholTry;
        return;
    end
end

error('%s is not positive definite even after numerical ridge. Use fewer metabolites or inspect the matrix.', matrixName);

end

% -------------------------------------------------------------------------
function p = ChiSquareUpperTail(T, df)

if ~isfinite(T) || ~isfinite(df) || df <= 0
    p = NaN;
    return;
end

try
    p = chi2cdf(T, df, 'upper');
catch
    p = 1 - chi2cdf(T, df);
end

end

% -------------------------------------------------------------------------
function [analysisData, lrtCfg] = ResolveCentralizedInputs(inputA, varargin)
if isstruct(inputA) && isfield(inputA, 'kind') && ...
        strcmp(string(inputA.kind), "MRSAnalysisData")
    analysisData = inputA;
    if isempty(varargin), lrtCfg = struct(); else, lrtCfg = varargin{1}; end
    return;
end

if isempty(varargin)
    error('Legacy interface requires covOutputs and deGraafOutputs.');
end
deGraafOutputs = varargin{1};
if numel(varargin) >= 2 && ~isempty(varargin{2})
    lrtCfg = varargin{2};
else
    lrtCfg = struct();
end

filterCfg = struct();
filterCfg.patientIDs = GetOption(lrtCfg, 'patientIDs', "all");
filterCfg.useSumPreferredFilter = GetOption(lrtCfg, 'excludeSumMetabolites', true);
filterCfg.sumMetabolites = GetOption(lrtCfg, 'sumMetabolites', ...
    ["GPC+PCh", "NAA+NAAG", "Cr+PCr", "Glu+Gln"]);
filterCfg.useCRLBMajorityFilter = GetOption(lrtCfg, 'useCRLBMajorityFilter', true);
filterCfg.crlbMajorityThreshold = GetOption(lrtCfg, 'crlbMajorityThreshold', 100);
filterCfg.ignoreZeros = GetOption(lrtCfg, 'ignoreZeros', true);
filterCfg.prepareTemporalCircularShift = false;
filterCfg.wishartViews.legacy = struct( ...
    'metabolites', GetOption(lrtCfg, 'metabolites', "all"), ...
    'minValidParts', GetOption(lrtCfg, 'minValidParts', 10));
[analysisData, ~] = ApplyAnalysisFilters(inputA, deGraafOutputs, filterCfg);
lrtCfg.filterView = "legacy";
end

% -------------------------------------------------------------------------
function [covOutputs, modelOutputs] = ExtractCanonicalOutputs(analysisData)
covOutputs = struct();
covOutputs.patientResultsByID = struct();
modelOutputs = struct();
modelOutputs.patientResultsByID = struct();

for p = 1:numel(analysisData.patientIDs)
    field = char(analysisData.patientIDs(p));
    entry = analysisData.patientDataByID.(field);
    covOutputs.patientResultsByID.(field).partTable = entry.partTable;
    modelOutputs.patientResultsByID.(field).meanAmplitudeCovTable = ...
        entry.modelCovarianceTable;
    modelOutputs.patientResultsByID.(field).meanAmplitudeCorrTable = ...
        entry.modelCorrelationTable;
end

end
