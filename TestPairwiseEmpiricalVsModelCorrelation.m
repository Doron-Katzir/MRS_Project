function pairOutputs = TestPairwiseEmpiricalVsModelCorrelation(inputA, varargin)
%TestPairwiseEmpiricalVsModelCorrelation Statistical layer for filtered data.
%
% Preferred interface:
%   outputs = TestPairwiseEmpiricalVsModelCorrelation(analysisData, statsCfg)
%
% The legacy (covOutputs, deGraafOutputs, pairCfg) interface is retained as
% an adapter; it delegates every inclusion decision to ApplyAnalysisFilters.

[analysisData, statsCfg] = ResolveInputs(inputA, varargin{:});
opts = ParseStatisticalOptions(statsCfg);
view = analysisData.pairwise;

fprintf('\nRunning pairwise empirical-vs-LCModel/deGraaf correlation tests...\n');
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

if opts.exportResults
    if ~exist(opts.outputDir, 'dir'), mkdir(opts.outputDir); end
    writetable(pairSummaryTable, fullfile(opts.outputDir, 'Pairwise_Empirical_vs_Model_Summary.csv'));
    writetable(patientPairTable, fullfile(opts.outputDir, 'Pairwise_Empirical_vs_Model_PatientLevel.csv'));
    fprintf('\nExported pairwise correlation test results to:\n  %s\n', opts.outputDir);
end
fprintf('\nPairwise empirical-vs-model correlation tests complete.\n');
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
