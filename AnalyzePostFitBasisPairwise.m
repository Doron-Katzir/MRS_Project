function outputs = AnalyzePostFitBasisPairwise( ...
        covOutputs, deGraafOutputs, analysisData, statsCfg, cfg, basisCfg)
%AnalyzePostFitBasisPairwise Parallel Division-1 post-fit basis experiment.
%
% Experiment A solves the independent-component covariance from B'*B.
% Experiment B appends the same fit's fitted baseline as one nuisance column,
% solves the full covariance, and then retains the component block. Both
% component covariances use the same reported-sum propagation and are averaged
% across parts before conversion to one patient correlation matrix. The
% already-filtered production pairwise view supplies the unchanged empirical
% correlations and eligibility policy; only rModel is replaced.

if nargin < 4 || isempty(statsCfg), statsCfg = struct(); end
if nargin < 5 || isempty(cfg), cfg = ProjectConfig(); end
if nargin < 6 || isempty(basisCfg), basisCfg = struct(); end
opts = ParseOptions(cfg, basisCfg);

ValidateTopLevelInputs(covOutputs, deGraafOutputs, analysisData, opts);
patientIDs = string(analysisData.pairwise.patientIDs(:));
nPatients = numel(patientIDs);
nReported = numel(opts.reportedNames);
nExpectedParts = numel(opts.expectedParts);

fprintf('\nPost-fit basis empirical-vs-model diagnostic\n');
fprintf('Patients: %d; Division: %d; expected parts per patient: %d\n', ...
    nPatients, opts.division, nExpectedParts);
fprintf('Experiment A: B''B; Experiment B: [B baseline]''[B baseline]\n');

patientResults = repmat(EmptyPatientResult(), nPatients, 1);
patientResultsByID = struct();
allPartDiagnostics = repmat(EmptyPartDiagnostic(), nPatients*nExpectedParts, 1);
diagnosticIndex = 0;

for pIdx = 1:nPatients
    patientID = patientIDs(pIdx);
    fprintf('  %s (%d/%d)\n', patientID, pIdx, nPatients);
    covPatient = GetPatientByID(covOutputs, patientID);
    modelPatient = GetPatientByID(deGraafOutputs, patientID);
    ValidatePatientScope(covPatient, modelPatient, patientID, opts.division);

    C0Stack = nan(nReported, nReported, nExpectedParts);
    C1Stack = nan(nReported, nReported, nExpectedParts);
    componentNamesByPart = cell(nExpectedParts, 1);

    for expectedIdx = 1:nExpectedParts
        partNumber = opts.expectedParts(expectedIdx);
        diagnosticIndex = diagnosticIndex + 1;
        d = EmptyPartDiagnostic();
        d.patientID = patientID;
        d.division = opts.division;
        d.part = partNumber;

        try
            partRow = FindPartRow(modelPatient.partInfo, partNumber);
            if isempty(partRow)
                d.status = "missing";
                d.errorMessage = "No parsed Division-1 fit is available for this expected part.";
                allPartDiagnostics(diagnosticIndex) = d;
                continue;
            end

            d.coordFile = string(partRow.coordFile);
            ValidateCoordFileScope(d.coordFile, opts.division, partNumber);
            [fitData, quantRows] = PostFitBasisCovariance.GetLoadedPartData( ...
                covPatient, d.coordFile);
            [B, componentNames, construction] = ...
                PostFitBasisCovariance.BuildIndependentBasis(fitData, quantRows);
            [M, ~] = PostFitBasisCovariance.BuildReportedTransformation( ...
                componentNames, opts.reportedNames, opts.independentCandidateNames);

            d.available = true;
            d.spectralPointCount = size(B, 1);
            d.nComponents = size(B, 2);
            d.componentNames = strjoin(componentNames, ", ");
            d.concentrationNormalized = true;
            componentNamesByPart{expectedIdx} = componentNames;

            [C0Component, h0] = PostFitBasisCovariance.SolveGeometry(B);
            d.rankNoBaseline = h0.rankH;
            d.conditionNumberNoBaseline = h0.conditionNumber;
            d.rcondNoBaseline = h0.rcondH;
            d.positiveDefiniteNoBaseline = h0.positiveDefinite;
            d.usableNoBaseline = h0.numericallyUsable;
            if h0.numericallyUsable
                C0 = M * C0Component * M.';
                C0Stack(:, :, expectedIdx) = (C0 + C0.') ./ 2;
            end

            baseline = PostFitBasisCovariance.GetBaselineColumn( ...
                fitData, size(B, 1));
            d.baselineFromSameFit = true;
            d.baselineNorm = norm(baseline);
            [C1Full, h1] = PostFitBasisCovariance.SolveGeometry([B, baseline]);
            d.rankWithBaseline = h1.rankH;
            d.conditionNumberWithBaseline = h1.conditionNumber;
            d.rcondWithBaseline = h1.rcondH;
            d.positiveDefiniteWithBaseline = h1.positiveDefinite;
            d.usableWithBaseline = h1.numericallyUsable;
            if h1.numericallyUsable
                nComponents = size(B, 2);
                C1Component = C1Full(1:nComponents, 1:nComponents);
                C1 = M * C1Component * M.';
                C1Stack(:, :, expectedIdx) = (C1 + C1.') ./ 2;
            end

            d.status = StatusFromUsability(d.usableNoBaseline, d.usableWithBaseline);
            if d.status ~= "ok"
                d.errorMessage = "At least one information matrix failed the shared numerical usability checks.";
            end
            d.nExcludedBasisTraces = height(construction.exclusionTable);
        catch ME
            d.status = "invalid";
            d.errorMessage = string(ME.message);
        end
        allPartDiagnostics(diagnosticIndex) = d;
    end

    C0Patient = MeanCovariance(C0Stack);
    C1Patient = MeanCovariance(C1Stack);
    R0Patient = PostFitBasisCovariance.CovarianceToCorrelation(C0Patient);
    R1Patient = PostFitBasisCovariance.CovarianceToCorrelation(C1Patient);
    patientPartDiagnostics = allPartDiagnostics( ...
        (diagnosticIndex-nExpectedParts+1):diagnosticIndex);

    result = EmptyPatientResult();
    result.patientID = patientID;
    result.division = opts.division;
    result.expectedParts = opts.expectedParts(:);
    result.reportedNames = opts.reportedNames;
    result.componentNamesByPart = componentNamesByPart;
    result.CNoBaselineByPart = C0Stack;
    result.CWithBaselineByPart = C1Stack;
    result.CNoBaselinePatient = C0Patient;
    result.CWithBaselinePatient = C1Patient;
    result.RNoBaselinePatient = R0Patient;
    result.RWithBaselinePatient = R1Patient;
    result.nPartsConsidered = nExpectedParts;
    result.nPartsAvailable = sum([patientPartDiagnostics.available]);
    result.nUsablePartsNoBaseline = sum([patientPartDiagnostics.usableNoBaseline]);
    result.nUsablePartsWithBaseline = sum([patientPartDiagnostics.usableWithBaseline]);
    result.missingParts = [patientPartDiagnostics( ...
        ~[patientPartDiagnostics.available]).part].';
    result.invalidPartsNoBaseline = [patientPartDiagnostics( ...
        [patientPartDiagnostics.available] & ...
        ~[patientPartDiagnostics.usableNoBaseline]).part].';
    result.invalidPartsWithBaseline = [patientPartDiagnostics( ...
        [patientPartDiagnostics.available] & ...
        ~[patientPartDiagnostics.usableWithBaseline]).part].';
    patientResults(pIdx) = result;
    patientResultsByID.(matlab.lang.makeValidName(char(patientID))) = result;
end

partDiagnosticsTable = struct2table(allPartDiagnostics);
patientSummaryTable = BuildPatientSummary(patientResults);
[C0Group, C1Group, R0Group, R1Group] = BuildGroupMatrices(patientResults);

noBaselineView = BuildBasisPairwiseView( ...
    analysisData.pairwise, patientResultsByID, opts.reportedNames, false);
withBaselineView = BuildBasisPairwiseView( ...
    analysisData.pairwise, patientResultsByID, opts.reportedNames, true);
AssertIdenticalEmpiricalViews(noBaselineView, withBaselineView, analysisData.pairwise);

pairStatsCfg = statsCfg;
pairStatsCfg.exportResults = false;
basisNoBaselinePairOutputs = RunExistingPairwiseTest( ...
    noBaselineView, analysisData, pairStatsCfg, "rBasisNoBaseline");
basisWithBaselinePairOutputs = RunExistingPairwiseTest( ...
    withBaselineView, analysisData, pairStatsCfg, "rBasisWithBaseline");

[comparisonTable, globalSummaryTable, directSummary] = CompareExperiments( ...
    basisNoBaselinePairOutputs, basisWithBaselinePairOutputs);

selectedNames = string(analysisData.pairwise.metabolites(:));
selectedIdx = FindNameIndices(opts.reportedNames, selectedNames);
group = struct();
group.reportedNames = opts.reportedNames;
group.CNoBaseline = C0Group;
group.CWithBaseline = C1Group;
group.RNoBaseline = R0Group;
group.RWithBaseline = R1Group;
group.selectedNames = selectedNames;
group.CNoBaselineSelected = C0Group(selectedIdx, selectedIdx);
group.CWithBaselineSelected = C1Group(selectedIdx, selectedIdx);
group.RNoBaselineSelected = R0Group(selectedIdx, selectedIdx);
group.RWithBaselineSelected = R1Group(selectedIdx, selectedIdx);
group.RDifferenceSelected = group.RWithBaselineSelected - group.RNoBaselineSelected;

outputFiles = struct();
if opts.exportResults
    outputFiles = ExportResults(basisNoBaselinePairOutputs, ...
        basisWithBaselinePairOutputs, comparisonTable, globalSummaryTable, ...
        partDiagnosticsTable, patientSummaryTable, group, opts);
end

figureFiles = strings(0, 1);
if opts.makeFigures
    figureFiles = MakeCohortFigures(group, comparisonTable, opts);
end

validation = BuildValidation(analysisData, patientResults, ...
    noBaselineView, withBaselineView, opts);
outputs = struct();
outputs.settings = opts;
outputs.patientIDs = patientIDs;
outputs.patientResults = patientResults;
outputs.patientResultsByID = patientResultsByID;
outputs.partDiagnosticsTable = partDiagnosticsTable;
outputs.patientSummaryTable = patientSummaryTable;
outputs.basisNoBaselinePairOutputs = basisNoBaselinePairOutputs;
outputs.basisWithBaselinePairOutputs = basisWithBaselinePairOutputs;
outputs.comparisonTable = comparisonTable;
outputs.globalSummaryTable = globalSummaryTable;
outputs.directSummary = directSummary;
outputs.group = group;
outputs.validation = validation;
outputs.outputFiles = outputFiles;
outputs.figureFiles = figureFiles;

PrintSummary(outputs);
end

% -------------------------------------------------------------------------
function opts = ParseOptions(cfg, basisCfg)
opts = struct();
opts.patientIDs = "all";
opts.division = 1;
opts.expectedPatientCount = 53;
opts.expectedParts = (1:36).';
opts.reportedNames = PostFitBasisCovariance.DefaultReportedNames();
opts.independentCandidateNames = ...
    PostFitBasisCovariance.DefaultIndependentCandidateNames();
opts.exportResults = true;
opts.makeFigures = true;
opts.figureVisible = "off";
opts.closeFigures = true;

rootDir = string(pwd);
if isstruct(cfg) && isfield(cfg, 'paths') && isfield(cfg.paths, 'rootDir')
    rootDir = string(cfg.paths.rootDir);
end
opts.noBaselineOutputDir = fullfile(rootDir, "PairwiseBasisModel_NoBaseline");
opts.withBaselineOutputDir = fullfile(rootDir, "PairwiseBasisModel_WithBaseline");
opts.comparisonOutputDir = fullfile(rootDir, "PairwiseBasisModel_Comparison");

provided = string(fieldnames(basisCfg));
allowed = string(fieldnames(opts));
for k = 1:numel(provided)
    if ~ismember(provided(k), allowed)
        error('Unknown post-fit basis pairwise option: %s', provided(k));
    end
    opts.(provided(k)) = basisCfg.(provided(k));
end
opts.patientIDs = string(opts.patientIDs(:));
opts.expectedParts = double(opts.expectedParts(:));
opts.reportedNames = string(opts.reportedNames(:));
opts.independentCandidateNames = string(opts.independentCandidateNames(:));
opts.figureVisible = string(opts.figureVisible);
if opts.division ~= 1
    error('This diagnostic is defined for Division 1 only.');
end
if ~isequal(opts.expectedParts, (1:36).')
    error('This diagnostic must consider exactly Division-1 parts 1:36.');
end
end

function ValidateTopLevelInputs(covOutputs, deGraafOutputs, analysisData, opts)
required = {'patientResultsByID'};
for k = 1:numel(required)
    if ~isfield(covOutputs, required{k}) || ~isfield(deGraafOutputs, required{k})
        error('Both loaded inputs must expose patientResultsByID.');
    end
end
if ~isstruct(analysisData) || ~isfield(analysisData, 'pairwise')
    error('analysisData must be the output of ApplyAnalysisFilters.');
end
if ~isfield(analysisData, 'division') || double(analysisData.division) ~= 1
    error('The centrally filtered analysisData must represent Division 1.');
end
ids = string(analysisData.pairwise.patientIDs(:));
if numel(ids) ~= opts.expectedPatientCount
    error('Expected exactly %d centrally selected patients; found %d.', ...
        opts.expectedPatientCount, numel(ids));
end
if ~(isscalar(opts.patientIDs) && strcmpi(opts.patientIDs, "all"))
    if ~isequal(ids, opts.patientIDs)
        error('Configured patientIDs do not match the central pairwise patient order.');
    end
end
end

function ValidatePatientScope(covPatient, modelPatient, patientID, division)
if double(covPatient.division) ~= division || double(modelPatient.division) ~= division
    error('Patient %s inputs do not both represent Division %d.', patientID, division);
end
end

function row = FindPartRow(partInfo, partNumber)
row = table();
if ~istable(partInfo) || isempty(partInfo), return; end
mask = double(partInfo.part) == partNumber & logical(partInfo.wasParsed);
idx = find(mask, 1, 'first');
if ~isempty(idx), row = partInfo(idx, :); end
end

function ValidateCoordFileScope(coordFile, division, partNumber)
[~, base, ext] = fileparts(char(coordFile));
token = regexp(string(base) + string(ext), ...
    'Division_(\d+)_(?:part_)?(\d+)\.basis\.coord$', 'tokens', 'once');
if isempty(token) || str2double(token{1}) ~= division || ...
        str2double(token{2}) ~= partNumber
    error('Coord file does not match Division %d part %d: %s', ...
        division, partNumber, coordFile);
end
end

function status = StatusFromUsability(noBaseline, withBaseline)
if noBaseline && withBaseline
    status = "ok";
elseif noBaseline
    status = "with_baseline_invalid";
elseif withBaseline
    status = "no_baseline_invalid";
else
    status = "both_invalid";
end
end

function C = MeanCovariance(stack)
C = mean(stack, 3, 'omitnan');
C = (C + C.') ./ 2;
end

function [C0, C1, R0, R1] = BuildGroupMatrices(patientResults)
nPatients = numel(patientResults);
n = size(patientResults(1).CNoBaselinePatient, 1);
stack0 = nan(n, n, nPatients);
stack1 = nan(n, n, nPatients);
for k = 1:nPatients
    stack0(:, :, k) = patientResults(k).CNoBaselinePatient;
    stack1(:, :, k) = patientResults(k).CWithBaselinePatient;
end
C0 = MeanCovariance(stack0);
C1 = MeanCovariance(stack1);
R0 = PostFitBasisCovariance.CovarianceToCorrelation(C0);
R1 = PostFitBasisCovariance.CovarianceToCorrelation(C1);
end

function view = BuildBasisPairwiseView(baseView, patientResultsByID, reportedNames, withBaseline)
view = baseView;
eligibility = baseView.patientEligibilityTable;
minParts = double(baseView.minValidParts);
for rowIdx = 1:height(eligibility)
    patientID = string(eligibility.patientID(rowIdx));
    metabA = string(eligibility.metaboliteA(rowIdx));
    metabB = string(eligibility.metaboliteB(rowIdx));
    patient = patientResultsByID.(matlab.lang.makeValidName(char(patientID)));
    if withBaseline
        R = patient.RWithBaselinePatient;
    else
        R = patient.RNoBaselinePatient;
    end
    rModel = GetNamedMatrixValue(R, reportedNames, metabA, metabB);
    eligibility.rModel(rowIdx) = rModel;

    if eligibility.nValidParts(rowIdx) < minParts
        eligibility.eligible(rowIdx) = false;
        eligibility.reason(rowIdx) = "insufficient valid parts";
    elseif ~isfinite(eligibility.rEmpirical(rowIdx))
        eligibility.eligible(rowIdx) = false;
        eligibility.reason(rowIdx) = "nonfinite empirical correlation";
    elseif ~isfinite(rModel)
        eligibility.eligible(rowIdx) = false;
        eligibility.reason(rowIdx) = "nonfinite model correlation";
    else
        eligibility.eligible(rowIdx) = true;
        eligibility.reason(rowIdx) = "";
    end
end
view.patientEligibilityTable = eligibility;
for pairIdx = 1:height(view.pairTable)
    mask = eligibility.metaboliteA == view.pairTable.metaboliteA(pairIdx) & ...
        eligibility.metaboliteB == view.pairTable.metaboliteB(pairIdx);
    nEligible = sum(eligibility.eligible(mask));
    view.pairTable.nEligiblePatients(pairIdx) = nEligible;
    view.pairTable.groupEligible(pairIdx) = ...
        nEligible >= view.minPatientsForGroupTest;
    if view.pairTable.groupEligible(pairIdx)
        view.pairTable.groupExclusionReason(pairIdx) = "";
    else
        view.pairTable.groupExclusionReason(pairIdx) = ...
            "fewer than " + string(view.minPatientsForGroupTest) + " eligible patients";
    end
end
end

function value = GetNamedMatrixValue(M, names, rowName, columnName)
r = find(names == rowName, 1, 'first');
c = find(names == columnName, 1, 'first');
if isempty(r) || isempty(c), value = NaN; else, value = M(r, c); end
end

function AssertIdenticalEmpiricalViews(view0, view1, productionView)
assert(isequal(view0.pairTable.metaboliteA, view1.pairTable.metaboliteA) && ...
    isequal(view0.pairTable.metaboliteB, view1.pairTable.metaboliteB), ...
    'A/B pair ordering changed.');
assert(isequaln(view0.patientEligibilityTable.rEmpirical, ...
    view1.patientEligibilityTable.rEmpirical), 'A/B empirical r values differ.');
assert(isequaln(view0.patientEligibilityTable.rEmpirical, ...
    productionView.patientEligibilityTable.rEmpirical), ...
    'Basis empirical r values differ from the production filtered view.');
assert(isequal(view0.patientEligibilityTable.nValidParts, ...
    productionView.patientEligibilityTable.nValidParts), ...
    'Basis empirical valid-part counts differ from production.');
end

function pairOutputs = RunExistingPairwiseTest(view, analysisData, statsCfg, modelName)
envelope = struct();
envelope.kind = "MRSAnalysisData";
envelope.pairwise = view;
envelope.crlbMajorityTable = analysisData.crlbMajorityTable;
envelope.filterReport = analysisData.filterReport;
pairOutputs = TestPairwiseEmpiricalVsModelCorrelation(envelope, statsCfg);
pairOutputs.modelCorrelationName = modelName;
pairOutputs.method = "Existing pairwise Fisher-z/t-test/BH-FDR logic with " + modelName;
end

function [comparison, globalSummary, direct] = CompareExperiments(noBaseline, withBaseline)
A = noBaseline.pairSummaryTable;
B = withBaseline.pairSummaryTable;
comparison = table();
for k = 1:height(A)
    idx = find(B.metaboliteA == A.metaboliteA(k) & ...
        B.metaboliteB == A.metaboliteB(k), 1, 'first');
    if isempty(idx), continue; end
    row = table();
    row.metaboliteA = A.metaboliteA(k);
    row.metaboliteB = A.metaboliteB(k);
    row.groupModelR_noBaseline = A.groupModelR(k);
    row.groupModelR_withBaseline = B.groupModelR(idx);
    row.meanDeltaZ_noBaseline = A.meanDeltaZ_empMinusModel(k);
    row.meanDeltaZ_withBaseline = B.meanDeltaZ_empMinusModel(idx);
    row.absMeanDeltaZ_noBaseline = abs(row.meanDeltaZ_noBaseline);
    row.absMeanDeltaZ_withBaseline = abs(row.meanDeltaZ_withBaseline);
    row.improvementAbsDeltaZ = row.absMeanDeltaZ_noBaseline - ...
        row.absMeanDeltaZ_withBaseline;
    row.nPatients_noBaseline = A.nPatientsUsed(k);
    row.nPatients_withBaseline = B.nPatientsUsed(idx);
    comparison = [comparison; row]; %#ok<AGROW>
end

common = isfinite(comparison.meanDeltaZ_noBaseline) & ...
    isfinite(comparison.meanDeltaZ_withBaseline);
a = comparison.meanDeltaZ_noBaseline(common);
b = comparison.meanDeltaZ_withBaseline(common);
globalSummary = [MakeGlobalRow("NoBaseline", a, A); ...
    MakeGlobalRow("WithBaseline", b, B)];
direct = struct();
direct.nCommonPairs = sum(common);
direct.fractionCloserWithBaseline = mean(abs(b) < abs(a));
direct.meanImprovementAbsDeltaZ = mean(abs(a) - abs(b));

validComparison = comparison(common, :);
[~, improveOrder] = sort(validComparison.improvementAbsDeltaZ, 'descend');
[~, worsenOrder] = sort(validComparison.improvementAbsDeltaZ, 'ascend');
nTop = min(10, height(validComparison));
direct.largestImprovements = validComparison(improveOrder(1:nTop), :);
direct.largestWorsenings = validComparison(worsenOrder(1:nTop), :);
end

function row = MakeGlobalRow(label, values, fullSummary)
ok = fullSummary.status == "ok" & isfinite(fullSummary.meanDeltaZ_empMinusModel);
row = table(label, numel(values), mean(values), median(values), std(values, 0), ...
    mean(abs(values)), sum(fullSummary.rejectH0_FDR & ok), ...
    'VariableNames', {'experiment','nPairs','meanPairMeanDeltaZ', ...
    'medianPairMeanDeltaZ','sdPairMeanDeltaZ','meanAbsPairMeanDeltaZ', ...
    'nFDRSignificantPairs'});
end

function files = ExportResults(noBaseline, withBaseline, comparison, globalSummary, ...
        partDiagnostics, patientSummary, group, opts)
EnsureDirectory(opts.noBaselineOutputDir);
EnsureDirectory(opts.withBaselineOutputDir);
EnsureDirectory(opts.comparisonOutputDir);
files = struct();
files.noBaselineSummaryCsv = string(fullfile(opts.noBaselineOutputDir, ...
    'Pairwise_BasisModel_NoBaseline_Summary.csv'));
files.noBaselinePatientCsv = string(fullfile(opts.noBaselineOutputDir, ...
    'Pairwise_BasisModel_NoBaseline_PatientLevel.csv'));
files.withBaselineSummaryCsv = string(fullfile(opts.withBaselineOutputDir, ...
    'Pairwise_BasisModel_WithBaseline_Summary.csv'));
files.withBaselinePatientCsv = string(fullfile(opts.withBaselineOutputDir, ...
    'Pairwise_BasisModel_WithBaseline_PatientLevel.csv'));
files.comparisonCsv = string(fullfile(opts.comparisonOutputDir, ...
    'Pairwise_BasisModel_AB_Comparison.csv'));
files.globalSummaryCsv = string(fullfile(opts.comparisonOutputDir, ...
    'Pairwise_BasisModel_AB_GlobalSummary.csv'));
files.partDiagnosticsCsv = string(fullfile(opts.comparisonOutputDir, ...
    'Pairwise_BasisModel_PartDiagnostics.csv'));
files.patientSummaryCsv = string(fullfile(opts.comparisonOutputDir, ...
    'Pairwise_BasisModel_PatientSummary.csv'));
files.groupNoBaselineCsv = string(fullfile(opts.comparisonOutputDir, ...
    'GroupModelCorrelation_NoBaseline.csv'));
files.groupWithBaselineCsv = string(fullfile(opts.comparisonOutputDir, ...
    'GroupModelCorrelation_WithBaseline.csv'));
files.groupDifferenceCsv = string(fullfile(opts.comparisonOutputDir, ...
    'GroupModelCorrelation_WithMinusNoBaseline.csv'));

writetable(noBaseline.pairSummaryTable, files.noBaselineSummaryCsv);
writetable(noBaseline.patientPairTable, files.noBaselinePatientCsv);
writetable(withBaseline.pairSummaryTable, files.withBaselineSummaryCsv);
writetable(withBaseline.patientPairTable, files.withBaselinePatientCsv);
writetable(comparison, files.comparisonCsv);
writetable(globalSummary, files.globalSummaryCsv);
writetable(partDiagnostics, files.partDiagnosticsCsv);
writetable(patientSummary, files.patientSummaryCsv);
writetable(MatrixTable(group.RNoBaselineSelected, group.selectedNames), ...
    files.groupNoBaselineCsv, 'WriteRowNames', true);
writetable(MatrixTable(group.RWithBaselineSelected, group.selectedNames), ...
    files.groupWithBaselineCsv, 'WriteRowNames', true);
writetable(MatrixTable(group.RDifferenceSelected, group.selectedNames), ...
    files.groupDifferenceCsv, 'WriteRowNames', true);
end

function figureFiles = MakeCohortFigures(group, comparison, opts)
EnsureDirectory(opts.comparisonOutputDir);
figureFiles = strings(5, 1);
figureFiles(1) = SaveMatrixFigure(group.RNoBaselineSelected, group.selectedNames, ...
    "Group post-fit basis model: no baseline", [-1, 1], ...
    fullfile(opts.comparisonOutputDir, 'Figure_1_GroupModel_NoBaseline.png'), opts);
figureFiles(2) = SaveMatrixFigure(group.RWithBaselineSelected, group.selectedNames, ...
    "Group post-fit basis model: with fitted baseline", [-1, 1], ...
    fullfile(opts.comparisonOutputDir, 'Figure_2_GroupModel_WithBaseline.png'), opts);
limit = max(abs(group.RDifferenceSelected), [], 'all', 'omitnan');
if ~isfinite(limit) || limit == 0, limit = 1; end
figureFiles(3) = SaveMatrixFigure(group.RDifferenceSelected, group.selectedNames, ...
    "Group model difference: with baseline - no baseline", [-limit, limit], ...
    fullfile(opts.comparisonOutputDir, 'Figure_3_GroupModel_Difference.png'), opts);

common = isfinite(comparison.meanDeltaZ_noBaseline) & ...
    isfinite(comparison.meanDeltaZ_withBaseline);
a = comparison.meanDeltaZ_noBaseline(common);
b = comparison.meanDeltaZ_withBaseline(common);
allValues = [a; b];
if isempty(allValues)
    edges = linspace(-1, 1, 21);
else
    lo = min(allValues); hi = max(allValues);
    if lo == hi, lo = lo - 0.5; hi = hi + 0.5; end
    edges = linspace(lo, hi, 21);
end
fig = figure('Visible', char(opts.figureVisible), 'Color', 'w');
histogram(a, edges, 'DisplayStyle', 'stairs', 'LineWidth', 2); hold on;
histogram(b, edges, 'DisplayStyle', 'stairs', 'LineWidth', 2);
xline(0, 'k:'); grid on;
xlabel('Pair-level mean \Delta z (empirical - model)'); ylabel('Number of pairs');
legend('No baseline', 'With baseline', 'Location', 'best');
title('Pairwise empirical-vs-model discrepancy');
figureFiles(4) = SaveFigure(fig, ...
    fullfile(opts.comparisonOutputDir, 'Figure_4_MeanDeltaZ_Histogram.png'), opts);

fig = figure('Visible', char(opts.figureVisible), 'Color', 'w');
x = comparison.groupModelR_noBaseline(common);
y = comparison.groupModelR_withBaseline(common);
scatter(x, y, 42, 'filled'); hold on; plot([-1, 1], [-1, 1], 'k--');
axis square; xlim([-1, 1]); ylim([-1, 1]); grid on;
xlabel('Group model r: no baseline'); ylabel('Group model r: with baseline');
title('Matched pair model correlations');
figureFiles(5) = SaveFigure(fig, ...
    fullfile(opts.comparisonOutputDir, 'Figure_5_GroupModelR_Scatter.png'), opts);
end

function file = SaveMatrixFigure(M, names, titleText, limits, filename, opts)
fig = figure('Visible', char(opts.figureVisible), 'Color', 'w');
imagesc(M); axis image; colorbar; clim(limits); colormap(BlueWhiteRedMap(256));
title(titleText, 'Interpreter', 'none');
xticks(1:numel(names)); xticklabels(names); xtickangle(45);
yticks(1:numel(names)); yticklabels(names);
set(gca, 'TickLabelInterpreter', 'none');
file = SaveFigure(fig, filename, opts);
end

function file = SaveFigure(fig, filename, opts)
exportgraphics(fig, filename, 'Resolution', 180);
file = string(filename);
if opts.closeFigures, close(fig); end
end

function map = BlueWhiteRedMap(n)
half = floor(n/2);
blueToWhite = [linspace(0,1,half).', linspace(0,1,half).', ones(half,1)];
redCount = n-half;
whiteToRed = [ones(redCount,1), linspace(1,0,redCount).', linspace(1,0,redCount).'];
map = [blueToWhite; whiteToRed];
end

function validation = BuildValidation(analysisData, patientResults, view0, view1, opts)
validation = struct();
validation.divisionOneOnly = opts.division == 1 && analysisData.division == 1;
validation.patientCount = numel(patientResults);
validation.expectedPatientCount = opts.expectedPatientCount;
validation.partsConsideredPerPatient = unique([patientResults.nPartsConsidered]);
validation.expectedParts = opts.expectedParts;
validation.sameBasisExtractionAndNormalization = true;
validation.onlyDesignDifferenceIsBaselineColumn = true;
validation.baselineAddedBeforeInversion = true;
validation.baselineBlockRemovedAfterInversion = true;
validation.sumsPropagatedAfterComponentCovariance = true;
validation.partCovarianceAveragedBeforeCorrelation = true;
validation.empiricalIdentical = isequaln( ...
    view0.patientEligibilityTable.rEmpirical, view1.patientEligibilityTable.rEmpirical);
validation.filterThresholdsIdentical = view0.minValidParts == view1.minValidParts && ...
    view0.minPatientsForGroupTest == view1.minPatientsForGroupTest;
validation.pairOrderIdentical = isequal(view0.pairTable(:, 1:2), view1.pairTable(:, 1:2));
validation.existingPairwiseStatisticsReused = true;
validation.existingBenjaminiHochbergReused = true;
validation.deGraafPrintCorrelationNotUsedAsBasisModel = true;
validation.noLCModelFitRunOrFileMutation = true;
end

function PrintSummary(outputs)
disp('Post-fit basis pairwise global summary:');
disp(outputs.globalSummaryTable);
fprintf('Common pairs: %d\n', outputs.directSummary.nCommonPairs);
fprintf('Fraction closer with baseline: %.6g\n', ...
    outputs.directSummary.fractionCloserWithBaseline);
fprintf('Mean improvement in |mean Delta-z|: %.6g\n', ...
    outputs.directSummary.meanImprovementAbsDeltaZ);
fprintf('Missing part rows: %d; invalid no-baseline rows: %d; invalid with-baseline rows: %d\n', ...
    sum(outputs.partDiagnosticsTable.status == "missing"), ...
    sum(outputs.partDiagnosticsTable.available & ~outputs.partDiagnosticsTable.usableNoBaseline), ...
    sum(outputs.partDiagnosticsTable.available & ~outputs.partDiagnosticsTable.usableWithBaseline));
end

function summary = BuildPatientSummary(patientResults)
summary = table(string({patientResults.patientID}).', ...
    [patientResults.nPartsConsidered].', [patientResults.nPartsAvailable].', ...
    [patientResults.nUsablePartsNoBaseline].', ...
    [patientResults.nUsablePartsWithBaseline].', ...
    string(arrayfun(@(p)JoinParts(p.missingParts), patientResults, 'UniformOutput', false)), ...
    string(arrayfun(@(p)JoinParts(p.invalidPartsNoBaseline), patientResults, 'UniformOutput', false)), ...
    string(arrayfun(@(p)JoinParts(p.invalidPartsWithBaseline), patientResults, 'UniformOutput', false)), ...
    'VariableNames', {'patientID','nPartsConsidered','nPartsAvailable', ...
    'nUsablePartsNoBaseline','nUsablePartsWithBaseline','missingParts', ...
    'invalidPartsNoBaseline','invalidPartsWithBaseline'});
end

function value = JoinParts(parts)
if isempty(parts), value = ""; else, value = strjoin(string(parts(:).'), ","); end
end

function patient = GetPatientByID(outputs, patientID)
field = matlab.lang.makeValidName(char(patientID));
if ~isfield(outputs.patientResultsByID, field)
    error('Patient %s is absent from a loaded input.', patientID);
end
patient = outputs.patientResultsByID.(field);
end

function idx = FindNameIndices(allNames, selectedNames)
[found, idx] = ismember(selectedNames, allNames);
if ~all(found)
    error('Selected metabolites absent from the configured reported panel: %s', ...
        strjoin(selectedNames(~found), ', '));
end
end

function T = MatrixTable(M, names)
T = array2table(M, 'VariableNames', matlab.lang.makeValidName(cellstr(names)), ...
    'RowNames', cellstr(names));
end

function EnsureDirectory(directory)
if ~exist(directory, 'dir'), mkdir(directory); end
end

function d = EmptyPartDiagnostic()
d = struct('patientID', "", 'division', NaN, 'part', NaN, 'status', "not_run", ...
    'available', false, 'coordFile', "", 'spectralPointCount', NaN, ...
    'nComponents', NaN, 'componentNames', "", 'concentrationNormalized', false, ...
    'nExcludedBasisTraces', NaN, 'baselineFromSameFit', false, 'baselineNorm', NaN, ...
    'rankNoBaseline', NaN, 'conditionNumberNoBaseline', NaN, 'rcondNoBaseline', NaN, ...
    'positiveDefiniteNoBaseline', false, 'usableNoBaseline', false, ...
    'rankWithBaseline', NaN, 'conditionNumberWithBaseline', NaN, ...
    'rcondWithBaseline', NaN, 'positiveDefiniteWithBaseline', false, ...
    'usableWithBaseline', false, 'errorMessage', "");
end

function p = EmptyPatientResult()
p = struct('patientID', "", 'division', NaN, 'expectedParts', zeros(0,1), ...
    'reportedNames', strings(0,1), 'componentNamesByPart', {cell(0,1)}, ...
    'CNoBaselineByPart', [], 'CWithBaselineByPart', [], ...
    'CNoBaselinePatient', [], 'CWithBaselinePatient', [], ...
    'RNoBaselinePatient', [], 'RWithBaselinePatient', [], ...
    'nPartsConsidered', 0, 'nPartsAvailable', 0, ...
    'nUsablePartsNoBaseline', 0, 'nUsablePartsWithBaseline', 0, ...
    'missingParts', zeros(0,1), 'invalidPartsNoBaseline', zeros(0,1), ...
    'invalidPartsWithBaseline', zeros(0,1));
end
