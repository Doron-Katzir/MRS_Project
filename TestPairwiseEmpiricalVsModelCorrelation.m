function pairOutputs = TestPairwiseEmpiricalVsModelCorrelation(inputA, varargin)
%TestPairwiseEmpiricalVsModelCorrelation Statistical layer for filtered data.
%
% Preferred interface:
%   outputs = TestPairwiseEmpiricalVsModelCorrelation(analysisData, statsCfg)
%
% The legacy (covOutputs, deGraafOutputs, pairCfg) interface is retained as
% an adapter; it delegates every inclusion decision to ApplyAnalysisFilters.
% statsCfg.empiricalSignalSource defaults to "original".  The alternate
% "baselinePCResidual" source consumes the already-aligned correlation table
% produced by PlotLCModelBaselineDiagnostics; it does not refit regressions.

[analysisData, statsCfg] = ResolveInputs(inputA, varargin{:});
opts = ParseStatisticalOptions(statsCfg);
view = analysisData.pairwise;
[view, empiricalSignalSource] = ResolveEmpiricalSignalView(view, statsCfg);

if empiricalSignalSource == "original"
    fprintf('\nRunning pairwise empirical-vs-LCModel/deGraaf correlation tests...\n');
else
    fprintf(['\nRunning baseline-PC-residual empirical-vs-LCModel/deGraaf ', ...
        'correlation tests...\n']);
end
fprintf('Number of patients available: %d\n', numel(view.patientIDs));
fprintf('Number of metabolites requested: %d\n', numel(view.metabolites));

pairSummaryTable = table();
patientPairTable = table();

for pairIdx = 1:height(view.pairTable)
    metabA = string(view.pairTable.metaboliteA(pairIdx));
    metabB = string(view.pairTable.metaboliteB(pairIdx));
    sourceRows = view.patientEligibilityTable.metaboliteA == metabA & ...
        view.patientEligibilityTable.metaboliteB == metabB;
    eligibility = view.patientEligibilityTable(sourceRows, :);
    thisPatientRows = table();

    for pIdx = 1:height(eligibility)
        source = eligibility(pIdx, :);
        row = MakeEmptyPatientPairRow(source.patientID, metabA, metabB);
        row.nValidParts = source.nValidParts;

        if source.eligible
            assert(isfinite(source.rEmpirical) && isfinite(source.rModel), ...
                'Central filtering invariant violated: eligible correlations must be finite.');
            rEmp = source.rEmpirical;
            rModel = source.rModel;
            zEmp = atanh(ClampCorrelation(rEmp));
            zModel = atanh(ClampCorrelation(rModel));

            row.status = "ok";
            row.rEmpirical = rEmp;
            row.rModel = rModel;
            row.zEmpirical = zEmp;
            row.zModel = zModel;
            row.deltaZ_empMinusModel = zEmp - zModel;
            row.deltaR_empMinusModel = rEmp - rModel;
            row.errorMessage = "";
        else
            row.status = "failed";
            row.errorMessage = FormatEligibilityReason(source, view.minValidParts);
        end
        thisPatientRows = [thisPatientRows; row]; %#ok<AGROW>
    end

    patientPairTable = [patientPairTable; thisPatientRows]; %#ok<AGROW>
    eligibleRows = eligibility.eligible;
    nPatientsUsed = sum(eligibleRows);
    summaryRow = MakeEmptyPairSummaryRow(metabA, metabB);
    summaryRow.nPatientsUsed = nPatientsUsed;

    if view.pairTable.groupEligible(pairIdx)
        assert(nPatientsUsed >= view.minPatientsForGroupTest, ...
            'Central filtering invariant violated: group eligibility is inconsistent.');
        d = thisPatientRows.deltaZ_empMinusModel(eligibleRows);
        zEmp = thisPatientRows.zEmpirical(eligibleRows);
        zModel = thisPatientRows.zModel(eligibleRows);
        rEmp = thisPatientRows.rEmpirical(eligibleRows);
        rModel = thisPatientRows.rModel(eligibleRows);
        nValidParts = thisPatientRows.nValidParts(eligibleRows);

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
        summaryRow.signConsistency = max(summaryRow.nPositiveDeltaZ, ...
            summaryRow.nNegativeDeltaZ) / nPatientsUsed;
        summaryRow.meanNvalidParts = mean(nValidParts);
        summaryRow.minNvalidParts = min(nValidParts);
        summaryRow.maxNvalidParts = max(nValidParts);
        summaryRow.errorMessage = "";

        fprintf('  %s vs %s: nPatients = %d, meanDeltaZ = %.4g, t = %.4g, p = %.4g\n', ...
            metabA, metabB, nPatientsUsed, summaryRow.meanDeltaZ_empMinusModel, tValue, pValue);
    else
        summaryRow.status = "failed";
        summaryRow.errorMessage = sprintf('Only %d valid patients; minimum requested is %d.', ...
            nPatientsUsed, view.minPatientsForGroupTest);
    end
    pairSummaryTable = [pairSummaryTable; summaryRow]; %#ok<AGROW>
end

pairSummaryTable.qValue_FDR = BenjaminiHochberg(pairSummaryTable.pValue);
pairSummaryTable.rejectH0_FDR = pairSummaryTable.qValue_FDR < opts.alpha;
sortP = pairSummaryTable.pValue;
sortP(~isfinite(sortP)) = Inf;
[~, order] = sort(sortP, 'ascend');
pairSummaryTable = pairSummaryTable(order, :);

pairOutputs = struct();
pairOutputs.method = "Pairwise Fisher-z empirical-vs-LCModel/deGraaf correlation t-test";
pairOutputs.nullHypothesis = "For each metabolite pair, mean over patients of atanh(r_empirical) - atanh(r_model) equals zero.";
pairOutputs.options = opts;
pairOutputs.patientIDs = view.patientIDs;
pairOutputs.metabolites = view.metabolites;
pairOutputs.crlbMajorityTable = analysisData.crlbMajorityTable;
pairOutputs.filterReport = analysisData.filterReport;
pairOutputs.pairSummaryTable = pairSummaryTable;
pairOutputs.patientPairTable = patientPairTable;
if empiricalSignalSource ~= "original"
    pairOutputs.empiricalSignalType = empiricalSignalSource;
end

if opts.exportResults
    if ~exist(opts.outputDir, 'dir'), mkdir(opts.outputDir); end
    if empiricalSignalSource == "original"
        summaryFile = 'Pairwise_Empirical_vs_Model_Summary.csv';
        patientFile = 'Pairwise_Empirical_vs_Model_PatientLevel.csv';
    else
        summaryFile = ...
            'Pairwise_Empirical_vs_Model_BaselineAdjusted_Summary.csv';
        patientFile = ...
            'Pairwise_Empirical_vs_Model_BaselineAdjusted_PatientLevel.csv';
    end
    writetable(pairSummaryTable, fullfile(opts.outputDir, summaryFile));
    writetable(patientPairTable, fullfile(opts.outputDir, patientFile));
    fprintf('\nExported pairwise correlation test results to:\n  %s\n', opts.outputDir);
end
if empiricalSignalSource == "original"
    fprintf('\nPairwise empirical-vs-model correlation tests complete.\n');
else
    fprintf('\nBaseline-adjusted pairwise empirical-vs-model tests complete.\n');
end
end

% -------------------------------------------------------------------------
function [analysisData, statsCfg] = ResolveInputs(inputA, varargin)
if IsAnalysisData(inputA)
    analysisData = inputA;
    if isempty(varargin), statsCfg = struct(); else, statsCfg = varargin{1}; end
    return;
end

if numel(varargin) < 1
    error('Legacy interface requires covOutputs and deGraafOutputs.');
end
deGraafOutputs = varargin{1};
if numel(varargin) >= 2 && ~isempty(varargin{2})
    legacyCfg = varargin{2};
else
    legacyCfg = struct();
end

filterCfg = struct();
filterCfg.patientIDs = GetOption(legacyCfg, 'patientIDs', "all");
filterCfg.metabolites = GetOption(legacyCfg, 'metabolites', "all");
filterCfg.useSumPreferredFilter = GetOption(legacyCfg, 'excludeSumMetabolites', true);
filterCfg.sumMetabolites = GetOption(legacyCfg, 'sumMetabolites', ...
    ["GPC+PCh", "NAA+NAAG", "Cr+PCr", "Glu+Gln"]);
filterCfg.useCRLBMajorityFilter = GetOption(legacyCfg, 'useCRLBMajorityFilter', true);
filterCfg.crlbMajorityThreshold = GetOption(legacyCfg, 'crlbMajorityThreshold', 100);
filterCfg.ignoreZeros = GetOption(legacyCfg, 'ignoreZeros', true);
filterCfg.pairwiseMinValidParts = GetOption(legacyCfg, 'minValidParts', 10);
filterCfg.pairwiseMinPatients = GetOption(legacyCfg, 'minPatientsForGroupTest', 3);
filterCfg.prepareTemporalCircularShift = false;
[analysisData, ~] = ApplyAnalysisFilters(inputA, deGraafOutputs, filterCfg);
statsCfg = legacyCfg;
end

function tf = IsAnalysisData(value)
tf = isstruct(value) && isfield(value, 'kind') && ...
    strcmp(string(value.kind), "MRSAnalysisData");
end

function opts = ParseStatisticalOptions(cfg)
opts = struct();
opts.alpha = double(GetOption(cfg, 'alpha', 0.05));
opts.exportResults = logical(GetOption(cfg, 'exportResults', false));
opts.outputDir = char(GetOption(cfg, 'outputDir', ...
    fullfile(pwd, 'PairwiseEmpiricalVsModelResults')));
end

function [view, signalSource] = ResolveEmpiricalSignalView(view, cfg)
signalSource = string(GetOption(cfg, 'empiricalSignalSource', "original"));
if ~isscalar(signalSource) || ismissing(signalSource) || ...
        strlength(signalSource) == 0
    error('empiricalSignalSource must be one nonempty string scalar.');
end

if strcmpi(signalSource, "original")
    signalSource = "original";
    return;
end
if ~strcmpi(signalSource, "baselinePCResidual")
    error('Unsupported empiricalSignalSource: %s.', signalSource);
end
signalSource = "baselinePCResidual";

if ~isstruct(cfg) || ~isfield(cfg, 'empiricalCorrelationTable') || ...
        ~istable(cfg.empiricalCorrelationTable)
    error(['baselinePCResidual requires statsCfg.empiricalCorrelationTable ', ...
        'from baselineDiagnostics.pairResidualizationTable.']);
end
view = SubstituteBaselineResidualCorrelations( ...
    view, cfg.empiricalCorrelationTable);
end

function view = SubstituteBaselineResidualCorrelations(view, residualTable)
requiredColumns = { ...
    'patientID','metaboliteA','metaboliteB','nValidPartsOriginal', ...
    'nValidPartsCommon','rEmpOriginal','rEmpResidual','rDeGraaf', ...
    'status','reason'};
missingColumns = requiredColumns(~ismember( ...
    requiredColumns, residualTable.Properties.VariableNames));
if ~isempty(missingColumns)
    error('The baseline residual table is missing column(s): %s.', ...
        strjoin(missingColumns, ', '));
end

eligibility = view.patientEligibilityTable;
sourceKeys = PairKeys(eligibility);
residualKeys = PairKeys(residualTable);
if numel(unique(residualKeys)) ~= numel(residualKeys)
    error('The baseline residual table contains duplicate patient/pair rows.');
end
[wasFound, locations] = ismember(sourceKeys, residualKeys);
if any(~wasFound)
    missingRow = find(~wasFound, 1, 'first');
    error('No baseline residual row was found for %s.', sourceKeys(missingRow));
end

originalEligibility = logical(eligibility.eligible);
for rowIndex = 1:height(eligibility)
    if ~originalEligibility(rowIndex)
        continue;
    end

    residualRow = residualTable(locations(rowIndex), :);
    if residualRow.nValidPartsOriginal ~= eligibility.nValidParts(rowIndex) || ...
            ~NearlyEqual(residualRow.rEmpOriginal, ...
            eligibility.rEmpirical(rowIndex), 1e-10) || ...
            ~NearlyEqual(residualRow.rDeGraaf, ...
            eligibility.rModel(rowIndex), 1e-12)
        error(['Baseline residual alignment check failed for patient %s, ', ...
            'pair %s/%s.'], eligibility.patientID(rowIndex), ...
            eligibility.metaboliteA(rowIndex), ...
            eligibility.metaboliteB(rowIndex));
    end

    eligibility.nValidParts(rowIndex) = residualRow.nValidPartsCommon;
    if string(residualRow.status) ~= "OK"
        eligibility.eligible(rowIndex) = false;
        eligibility.rEmpirical(rowIndex) = NaN;
        eligibility.reason(rowIndex) = ...
            "baseline-PC residual: " + string(residualRow.reason);
        continue;
    end
    if residualRow.nValidPartsCommon < view.minValidParts
        eligibility.eligible(rowIndex) = false;
        eligibility.rEmpirical(rowIndex) = NaN;
        eligibility.reason(rowIndex) = "insufficient valid parts";
        continue;
    end
    if ~isfinite(residualRow.rEmpResidual)
        eligibility.eligible(rowIndex) = false;
        eligibility.rEmpirical(rowIndex) = NaN;
        eligibility.reason(rowIndex) = ...
            "nonfinite baseline-PC residual correlation";
        continue;
    end

    eligibility.rEmpirical(rowIndex) = residualRow.rEmpResidual;
    eligibility.reason(rowIndex) = "";
end

view.patientEligibilityTable = eligibility;
for pairIndex = 1:height(view.pairTable)
    pairRows = eligibility.metaboliteA == ...
        string(view.pairTable.metaboliteA(pairIndex)) & ...
        eligibility.metaboliteB == ...
        string(view.pairTable.metaboliteB(pairIndex));
    nEligible = sum(eligibility.eligible(pairRows));
    view.pairTable.nEligiblePatients(pairIndex) = nEligible;
    view.pairTable.groupEligible(pairIndex) = ...
        nEligible >= view.minPatientsForGroupTest;
    if view.pairTable.groupEligible(pairIndex)
        view.pairTable.groupExclusionReason(pairIndex) = "";
    else
        view.pairTable.groupExclusionReason(pairIndex) = ...
            "fewer than " + string(view.minPatientsForGroupTest) + ...
            " eligible patients";
    end
end
end

function keys = PairKeys(T)
separator = string(char(31));
keys = string(T.patientID) + separator + string(T.metaboliteA) + ...
    separator + string(T.metaboliteB);
end

function tf = NearlyEqual(a, b, tolerance)
tf = isfinite(a) && isfinite(b) && abs(double(a) - double(b)) <= tolerance;
end

function value = GetOption(s, fieldName, defaultValue)
if isstruct(s) && isfield(s, fieldName), value = s.(fieldName); else, value = defaultValue; end
end

function message = FormatEligibilityReason(source, minValidParts)
switch string(source.reason)
    case "insufficient valid parts"
        message = string(sprintf('Only %d valid repeated parts; minimum requested is %d.', ...
            source.nValidParts, minValidParts));
    case "nonfinite empirical correlation"
        message = "Empirical correlation is non-finite.";
    case "nonfinite model correlation"
        message = "Model correlation is non-finite.";
    otherwise
        message = string(source.reason);
end
end

function row = MakeEmptyPatientPairRow(patientID, metabA, metabB)
row = table(string(patientID), string(metabA), string(metabB), "not_run", ...
    NaN, NaN, NaN, NaN, NaN, NaN, NaN, "", 'VariableNames', ...
    {'patientID','metaboliteA','metaboliteB','status','nValidParts', ...
    'rEmpirical','rModel','zEmpirical','zModel','deltaZ_empMinusModel', ...
    'deltaR_empMinusModel','errorMessage'});
end

function row = MakeEmptyPairSummaryRow(metabA, metabB)
row = table();
row.metaboliteA = string(metabA); row.metaboliteB = string(metabB); row.status = "not_run";
row.nPatientsUsed = NaN; row.df = NaN; row.meanDeltaZ_empMinusModel = NaN;
row.sdDeltaZ_empMinusModel = NaN; row.tValue = NaN; row.pValue = NaN; row.rejectH0 = false;
row.groupEmpiricalR = NaN; row.groupModelR = NaN; row.meanDeltaR_empMinusModel = NaN;
row.meanEmpiricalR = NaN; row.meanModelR = NaN; row.meanAbsEmpiricalR = NaN;
row.meanAbsModelR = NaN; row.meanDeltaAbsR_empMinusModel = NaN;
row.nPositiveDeltaZ = NaN; row.nNegativeDeltaZ = NaN; row.signConsistency = NaN;
row.meanNvalidParts = NaN; row.minNvalidParts = NaN; row.maxNvalidParts = NaN;
row.errorMessage = "";
end

function r = ClampCorrelation(r)
r = max(min(r, 0.999999), -0.999999);
end

function [tValue, pValue, df] = OneSampleTTestAgainstZero(x)
x = double(x(:));
x = x(isfinite(x));
n = numel(x); df = n - 1;
if n < 3, tValue = NaN; pValue = NaN; return; end
sx = std(x, 0);
if sx <= 0 || ~isfinite(sx), tValue = NaN; pValue = NaN; return; end
tValue = mean(x) / (sx / sqrt(n));
try
    pValue = 2 * tcdf(-abs(tValue), df);
catch
    pValue = NaN;
    warning('tcdf is unavailable. pValue was set to NaN.');
end
end

function q = BenjaminiHochberg(p)
p = p(:); q = nan(size(p)); valid = isfinite(p); pValid = p(valid);
if isempty(pValid), return; end
[pSorted, sortIdx] = sort(pValid, 'ascend'); m = numel(pSorted);
qSorted = pSorted .* m ./ (1:m)';
for k = m-1:-1:1, qSorted(k) = min(qSorted(k), qSorted(k+1)); end
qSorted = min(qSorted, 1); qValid = nan(size(pValid));
qValid(sortIdx) = qSorted; q(valid) = qValid;
end
