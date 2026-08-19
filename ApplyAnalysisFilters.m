function [analysisData, filterReport] = ApplyAnalysisFilters(empiricalOutputs, modelOutputs, filterCfg)
%ApplyAnalysisFilters Central scientific filtering for MRS statistical tests.
%
% Readers remain responsible for malformed or missing files. This function
% owns analysis choices: patient/metabolite selection, sum preference, CRLB
% majority, zero/nonfinite handling, and patient-by-pair eligibility.

if nargin < 3 || isempty(filterCfg)
    filterCfg = struct();
end

opts = ParseFilterOptions(filterCfg);
[patientIDs, patientReport] = ResolvePatients(empiricalOutputs, modelOutputs, opts.patientIDs);
startingMetabolites = ResolveStartingMetabolites(empiricalOutputs, opts.metabolites);

[patientDataByID, observationReport] = BuildCanonicalPatientData( ...
    empiricalOutputs, modelOutputs, patientIDs, opts.ignoreZeros);

crlbMajorityTable = BuildGlobalCRLBMajorityTable( ...
    modelOutputs, patientIDs, opts.crlbMajorityThreshold);

[selectedMetabolites, metaboliteReport] = ApplyMetabolitePolicy( ...
    startingMetabolites, "all", patientDataByID, patientIDs, ...
    crlbMajorityTable, opts);

[pairTable, patientPairEligibility] = BuildPairwiseEligibility( ...
    patientDataByID, patientIDs, selectedMetabolites, opts.pairwiseMinValidParts);

pairTable.groupEligible = pairTable.nEligiblePatients >= opts.pairwiseMinPatients;
pairTable.groupExclusionReason = strings(height(pairTable), 1);
pairTable.groupExclusionReason(~pairTable.groupEligible) = ...
    "fewer than " + string(opts.pairwiseMinPatients) + " eligible patients";

temporalMetabolites = startingMetabolites;
if opts.temporalUseGlobalMetabolites
    temporalMetabolites = selectedMetabolites;
end
[temporalPairTable, temporalPatientEligibility] = BuildEmpiricalPairEligibility( ...
    patientDataByID, patientIDs, temporalMetabolites, opts.temporalMinValidParts, ...
    opts.prepareTemporalCircularShift);
temporalPairTable.groupEligible = ...
    temporalPairTable.nEligiblePatients >= opts.temporalMinPatients;
temporalPairTable.shiftGroupEligible = ...
    temporalPairTable.nShiftEligiblePatients >= opts.temporalMinPatients;
temporalCRLBQualityTable = BuildTemporalCRLBQualityTable( ...
    patientDataByID, patientIDs, temporalMetabolites, opts.division, ...
    opts.temporalCRLBThreshold, opts.temporalRequiredGoodFraction);

wishartViews = BuildWishartViews(startingMetabolites, patientDataByID, ...
    patientIDs, crlbMajorityTable, opts);
canonicalEmpiricalGroup = BuildCanonicalEmpiricalGroup( ...
    patientDataByID, patientIDs, startingMetabolites, opts.empiricalSummaryMinValidPairs);

analysisData = struct();
analysisData.kind = "MRSAnalysisData";
analysisData.schemaVersion = 1;
analysisData.filterConfig = opts;
analysisData.division = opts.division;
analysisData.patientIDs = patientIDs;
analysisData.startingMetabolites = startingMetabolites;
analysisData.metabolites = selectedMetabolites;
analysisData.patientDataByID = patientDataByID;
analysisData.empiricalGroup = GetOptionalField(empiricalOutputs, 'group', struct());
analysisData.empiricalGroup.meanCorrMatrix = canonicalEmpiricalGroup.meanCorrMatrix;
analysisData.empiricalGroup.meanAbsCorrMatrix = canonicalEmpiricalGroup.meanAbsCorrMatrix;
analysisData.empiricalGroup.meanCorrTable = canonicalEmpiricalGroup.meanCorrTable;
analysisData.empiricalGroup.meanAbsCorrTable = canonicalEmpiricalGroup.meanAbsCorrTable;
analysisData.modelGroup = GetOptionalField(modelOutputs, 'group', struct());

analysisData.pairwise = struct();
analysisData.pairwise.patientIDs = patientIDs;
analysisData.pairwise.metabolites = selectedMetabolites;
analysisData.pairwise.pairTable = pairTable;
analysisData.pairwise.patientEligibilityTable = patientPairEligibility;
analysisData.pairwise.minValidParts = opts.pairwiseMinValidParts;
analysisData.pairwise.minPatientsForGroupTest = opts.pairwiseMinPatients;

analysisData.temporal = struct();
analysisData.temporal.patientIDs = patientIDs;
analysisData.temporal.metabolites = temporalMetabolites;
analysisData.temporal.pairTable = temporalPairTable;
analysisData.temporal.patientEligibilityTable = temporalPatientEligibility;
analysisData.temporal.minValidParts = opts.temporalMinValidParts;
analysisData.temporal.minPatients = opts.temporalMinPatients;
analysisData.temporal.crlbThreshold = opts.temporalCRLBThreshold;
analysisData.temporal.requiredGoodFraction = opts.temporalRequiredGoodFraction;
analysisData.temporal.crlbQualityTable = temporalCRLBQualityTable;
analysisData.temporal.circularShiftPrepared = opts.prepareTemporalCircularShift;
% Focused group summaries used only for the temporal result's descriptive
% empirical-vs-model columns. Raw patient/project state stays outside this view.
analysisData.temporal.empiricalMeanAbsCorrTable = ...
    analysisData.empiricalGroup.meanAbsCorrTable;
analysisData.temporal.modelMeanAmplitudeCorrTable = GetOptionalField( ...
    analysisData.modelGroup, 'meanAmplitudeCorrTable', table());
analysisData.temporal.modelMeanAbsAmplitudeCorrTable = GetOptionalField( ...
    analysisData.modelGroup, 'meanAbsAmplitudeCorrTable', table());

analysisData.wishart = struct();
analysisData.wishart.patientIDs = patientIDs;
analysisData.wishart.views = wishartViews;

analysisData.crlbMajorityTable = crlbMajorityTable;

filterReport = struct();
filterReport.config = opts;
filterReport.patientSummary = patientReport;
filterReport.startingPatientCount = patientReport.nEmpiricalPatients;
filterReport.commonPatientCount = numel(patientIDs);
filterReport.selectedPatientIDs = patientIDs;
filterReport.startingMetaboliteCount = numel(startingMetabolites);
filterReport.startingMetabolites = startingMetabolites;
filterReport.metaboliteDecisions = metaboliteReport;
filterReport.metabolitesRemovedBySum = metaboliteReport.metabolite( ...
    metaboliteReport.reason == "component of available preferred sum");
filterReport.metabolitesRemovedByCRLB = metaboliteReport.metabolite( ...
    metaboliteReport.reason == "failed global CRLB-majority criterion");
filterReport.finalMetaboliteCount = numel(selectedMetabolites);
filterReport.finalMetabolites = selectedMetabolites;
filterReport.crlbMajorityTable = crlbMajorityTable;
filterReport.observationSummary = observationReport;
filterReport.nZeroObservationsExcluded = sum(observationReport.nZerosConverted);
filterReport.nNonfiniteObservationsExcluded = sum(observationReport.nNonfiniteInput);
filterReport.patientPairEligibility = patientPairEligibility;
filterReport.pairSummary = pairTable;
filterReport.nPatientPairInsufficientValidParts = sum( ...
    patientPairEligibility.reason == "insufficient valid parts");

analysisData.filterReport = filterReport;
end

% -------------------------------------------------------------------------
function opts = ParseFilterOptions(cfg)
opts = struct();
opts.patientIDs = GetOption(cfg, 'patientIDs', "all");
opts.division = double(GetOption(cfg, 'division', 1));
opts.metabolites = string(GetOption(cfg, 'metabolites', "all"));
opts.useSumPreferredFilter = logical(GetOption(cfg, 'useSumPreferredFilter', true));
opts.sumMetabolites = string(GetOption(cfg, 'sumMetabolites', ...
    ["GPC+PCh", "NAA+NAAG", "Cr+PCr", "Glu+Gln"]));
opts.useCRLBMajorityFilter = logical(GetOption(cfg, 'useCRLBMajorityFilter', true));
opts.crlbMajorityThreshold = double(GetOption(cfg, 'crlbMajorityThreshold', 100));
opts.ignoreZeros = logical(GetOption(cfg, 'ignoreZeros', true));
opts.pairwiseMinValidParts = double(GetOption(cfg, 'pairwiseMinValidParts', 10));
opts.pairwiseMinPatients = double(GetOption(cfg, 'pairwiseMinPatients', 3));
opts.temporalMinValidParts = double(GetOption(cfg, 'temporalMinValidParts', 8));
opts.temporalMinPatients = double(GetOption(cfg, 'temporalMinPatients', 3));
opts.temporalUseGlobalMetabolites = logical(GetOption(cfg, 'temporalUseGlobalMetabolites', false));
opts.temporalCRLBThreshold = double(GetOption(cfg, 'temporalCRLBThreshold', 100));
opts.temporalRequiredGoodFraction = double(GetOption(cfg, 'temporalRequiredGoodFraction', 0.90));
opts.prepareTemporalCircularShift = logical(GetOption(cfg, 'prepareTemporalCircularShift', true));
opts.wishartMinValidParts = double(GetOption(cfg, 'wishartMinValidParts', 30));
opts.wishartViews = GetOption(cfg, 'wishartViews', struct());
opts.empiricalSummaryMinValidPairs = double(GetOption(cfg, 'empiricalSummaryMinValidPairs', 3));

if opts.division ~= 1
    error(['Current upstream outputs contain one selected division. ', ...
        'Behavior-preserving filtering requires division = 1.']);
end
end

% -------------------------------------------------------------------------
function [patientIDs, report] = ResolvePatients(empiricalOutputs, modelOutputs, requested)
if ~isfield(empiricalOutputs, 'patientResultsByID') || ...
        ~isfield(modelOutputs, 'patientResultsByID')
    error('Both inputs must contain patientResultsByID.');
end

empiricalIDs = string(fieldnames(empiricalOutputs.patientResultsByID));
modelIDs = string(fieldnames(modelOutputs.patientResultsByID));
commonIDs = intersect(empiricalIDs, modelIDs, 'stable');

if isnumeric(requested)
    requested = compose("P%02d", requested(:));
end
requested = string(requested(:));

if isscalar(requested) && strcmpi(requested, "all")
    patientIDs = commonIDs;
else
    requestedFields = strings(numel(requested), 1);
    for k = 1:numel(requested)
        requestedFields(k) = string(matlab.lang.makeValidName(char(requested(k))));
    end
    patientIDs = requestedFields(ismember(requestedFields, commonIDs));
end

if isempty(patientIDs)
    error('No common selected patients remain after filtering.');
end

report = struct();
report.nEmpiricalPatients = numel(empiricalIDs);
report.nModelPatients = numel(modelIDs);
report.nCommonBeforeSelection = numel(commonIDs);
report.nSelectedCommonPatients = numel(patientIDs);
report.empiricalOnlyPatientIDs = setdiff(empiricalIDs, modelIDs, 'stable');
report.modelOnlyPatientIDs = setdiff(modelIDs, empiricalIDs, 'stable');
report.commonPatientIDs = commonIDs;
end

% -------------------------------------------------------------------------
function metabs = ResolveStartingMetabolites(empiricalOutputs, requested)
if isfield(empiricalOutputs, 'metabList')
    available = string(empiricalOutputs.metabList(:));
elseif isfield(empiricalOutputs, 'group') && ...
        isfield(empiricalOutputs.group, 'meanCorrTable')
    T = empiricalOutputs.group.meanCorrTable;
    if ~isempty(T.Properties.RowNames)
        available = string(T.Properties.RowNames);
    else
        available = string(T.Properties.VariableNames);
    end
else
    error('Could not determine the empirical metabolite list.');
end

requested = string(requested(:));
if isscalar(requested) && strcmpi(requested, "all")
    metabs = available;
else
    metabs = KeepMatchingNames(available, requested);
end
metabs = metabs(:);
end

% -------------------------------------------------------------------------
function [dataByID, report] = BuildCanonicalPatientData(empiricalOutputs, modelOutputs, patientIDs, ignoreZeros)
dataByID = struct();
report = table();

for pIdx = 1:numel(patientIDs)
    patientID = patientIDs(pIdx);
    field = char(matlab.lang.makeValidName(char(patientID)));
    empiricalPatient = empiricalOutputs.patientResultsByID.(field);
    modelPatient = modelOutputs.patientResultsByID.(field);

    if ~isfield(empiricalPatient, 'partTable')
        error('Missing empirical partTable for patient %s.', patientID);
    end

    original = empiricalPatient.partTable;
    filtered = original;
    nZeros = 0;
    nNonfinite = 0;
    vars = string(filtered.Properties.VariableNames);
    for vIdx = 1:numel(vars)
        name = char(vars(vIdx));
        if strcmpi(name, 'part') || ~isnumeric(filtered.(name))
            continue;
        end
        values = double(filtered.(name));
        nNonfinite = nNonfinite + sum(~isfinite(values));
        if ignoreZeros
            zeroMask = values == 0;
            nZeros = nZeros + sum(zeroMask);
            values(zeroMask) = NaN;
        end
        values(~isfinite(values)) = NaN;
        filtered.(name) = values;
    end

    entry = struct();
    entry.patientID = patientID;
    entry.partTable = filtered;
    entry.modelCorrelationTable = GetOptionalField(modelPatient, 'meanAmplitudeCorrTable', table());
    entry.modelCovarianceTable = GetOptionalField(modelPatient, 'meanAmplitudeCovTable', table());
    entry.partCRLB = GetOptionalField(modelPatient, 'partCRLB', []);
    entry.metabList = string(GetOptionalField(modelPatient, 'metabList', strings(0, 1)));
    entry.empiricalCoordTable = GetOptionalField(empiricalPatient, 'coordTable', table());
    dataByID.(field) = entry;

    report = [report; table(patientID, height(original), nZeros, nNonfinite, ...
        'VariableNames', {'patientID','nParts','nZerosConverted','nNonfiniteInput'})]; %#ok<AGROW>
end
end

% -------------------------------------------------------------------------
function [metabs, decisions] = ApplyMetabolitePolicy(startingMetabs, requested, dataByID, patientIDs, crlbTable, opts)
startingMetabs = string(startingMetabs(:));
requested = string(requested(:));

if isscalar(requested) && strcmpi(requested, "all")
    metabs = startingMetabs;
else
    metabs = KeepMatchingNames(startingMetabs, requested);
end

decisions = table(startingMetabs, repmat("removed", numel(startingMetabs), 1), ...
    repmat("not requested", numel(startingMetabs), 1), ...
    'VariableNames', {'metabolite','status','reason'});
for k = 1:numel(metabs)
    idx = FindMatchingNameIndex(decisions.metabolite, metabs(k));
    decisions.status(idx) = "candidate";
    decisions.reason(idx) = "";
end

if opts.useSumPreferredFilter
    for sIdx = 1:numel(opts.sumMetabolites)
        sumName = opts.sumMetabolites(sIdx);
        if isempty(FindMatchingNameIndex(metabs, sumName))
            continue;
        end
        components = string(split(sumName, '+'));
        remove = false(numel(metabs), 1);
        for mIdx = 1:numel(metabs)
            if any(arrayfun(@(c) NamesMatch(metabs(mIdx), components(c)), 1:numel(components)))
                remove(mIdx) = true;
                dIdx = FindMatchingNameIndex(decisions.metabolite, metabs(mIdx));
                decisions.status(dIdx) = "removed";
                decisions.reason(dIdx) = "component of available preferred sum";
            end
        end
        metabs = metabs(~remove);
    end
end

if opts.useCRLBMajorityFilter
    allowed = crlbTable.metabolite(crlbTable.keepByCRLB);
    for mIdx = 1:numel(metabs)
        if isempty(FindMatchingNameIndex(allowed, metabs(mIdx)))
            dIdx = FindMatchingNameIndex(decisions.metabolite, metabs(mIdx));
            decisions.status(dIdx) = "removed";
            decisions.reason(dIdx) = "failed global CRLB-majority criterion";
        end
    end
    metabs = KeepMatchingNames(metabs, allowed);
end

available = false(numel(metabs), 1);
for mIdx = 1:numel(metabs)
    for pIdx = 1:numel(patientIDs)
        entry = dataByID.(char(patientIDs(pIdx)));
        if HasMatchingTableVariable(entry.partTable, metabs(mIdx)) && ...
                HasMatchingSquareLabel(entry.modelCorrelationTable, metabs(mIdx))
            available(mIdx) = true;
            break;
        end
    end
    if ~available(mIdx)
        dIdx = FindMatchingNameIndex(decisions.metabolite, metabs(mIdx));
        decisions.status(dIdx) = "removed";
        decisions.reason(dIdx) = "not available in empirical and model tables";
    end
end
metabs = metabs(available);

for mIdx = 1:numel(metabs)
    dIdx = FindMatchingNameIndex(decisions.metabolite, metabs(mIdx));
    decisions.status(dIdx) = "kept";
    decisions.reason(dIdx) = "";
end
end

% -------------------------------------------------------------------------
function views = BuildWishartViews(startingMetabs, dataByID, patientIDs, crlbTable, opts)
views = struct();
viewCfg = opts.wishartViews;
if isempty(fieldnames(viewCfg))
    viewCfg.default = struct('metabolites', "all");
end

names = string(fieldnames(viewCfg));
for k = 1:numel(names)
    name = char(names(k));
    cfg = viewCfg.(name);
    requested = GetOption(cfg, 'metabolites', "all");
    [metabs, decisions] = ApplyMetabolitePolicy(startingMetabs, requested, ...
        dataByID, patientIDs, crlbTable, opts);
    view = struct();
    view.patientIDs = patientIDs;
    view.metabolites = metabs;
    view.metaboliteDecisions = decisions;
    view.minValidParts = double(GetOption(cfg, 'minValidParts', opts.wishartMinValidParts));
    views.(name) = view;
end
end

% -------------------------------------------------------------------------
function group = BuildCanonicalEmpiricalGroup(dataByID, patientIDs, metabs, minValidPairs)
metabs = string(metabs(:));
stack = nan(numel(metabs), numel(metabs), numel(patientIDs));
for p = 1:numel(patientIDs)
    entry = dataByID.(char(patientIDs(p)));
    for a = 1:numel(metabs)
        x = GetTableColumn(entry.partTable, metabs(a));
        for b = 1:numel(metabs)
            y = GetTableColumn(entry.partTable, metabs(b));
            valid = isfinite(x) & isfinite(y);
            if sum(valid) >= minValidPairs
                stack(a, b, p) = PearsonR(x(valid), y(valid));
            end
        end
    end
end

meanCorr = nan(numel(metabs));
for a = 1:numel(metabs)
    for b = 1:numel(metabs)
        values = squeeze(stack(a, b, :));
        values = values(isfinite(values));
        if isempty(values), continue; end
        if a == b
            meanCorr(a, b) = 1;
        else
            values = max(min(values, 0.999999), -0.999999);
            meanCorr(a, b) = tanh(mean(atanh(values)));
        end
    end
end
meanAbs = mean(abs(stack), 3, 'omitnan');
group = struct();
group.meanCorrMatrix = meanCorr;
group.meanAbsCorrMatrix = meanAbs;
group.meanCorrTable = MatrixToMetaboliteTable(meanCorr, metabs);
group.meanAbsCorrTable = MatrixToMetaboliteTable(meanAbs, metabs);
end

function T = MatrixToMetaboliteTable(M, metabs)
T = array2table(M, 'VariableNames', ...
    matlab.lang.makeValidName(cellstr(metabs)), 'RowNames', cellstr(metabs));
end

% -------------------------------------------------------------------------
function [pairTable, patientTable] = BuildPairwiseEligibility(dataByID, patientIDs, metabs, minValidParts)
pairTable = table();
patientTable = table();
for a = 1:numel(metabs)
    for b = a+1:numel(metabs)
        metabA = metabs(a);
        metabB = metabs(b);
        firstRow = height(patientTable) + 1;
        for p = 1:numel(patientIDs)
            patientID = patientIDs(p);
            entry = dataByID.(char(patientID));
            [rEmp, rModel, nValid, eligible, reason] = EvaluatePair( ...
                entry, metabA, metabB, minValidParts, true);
            patientTable = [patientTable; table(patientID, metabA, metabB, ...
                nValid, eligible, rEmp, rModel, reason, ...
                'VariableNames', {'patientID','metaboliteA','metaboliteB', ...
                'nValidParts','eligible','rEmpirical','rModel','reason'})]; %#ok<AGROW>
        end
        rows = firstRow:height(patientTable);
        pairTable = [pairTable; table(metabA, metabB, sum(patientTable.eligible(rows)), ...
            'VariableNames', {'metaboliteA','metaboliteB','nEligiblePatients'})]; %#ok<AGROW>
    end
end
end

% -------------------------------------------------------------------------
function [pairTable, patientTable] = BuildEmpiricalPairEligibility(dataByID, patientIDs, metabs, minValidParts, prepareShifts)
pairTable = table();
patientTable = table();
for a = 1:numel(metabs)
    for b = a+1:numel(metabs)
        metabA = metabs(a);
        metabB = metabs(b);
        firstRow = height(patientTable) + 1;
        for p = 1:numel(patientIDs)
            patientID = patientIDs(p);
            entry = dataByID.(char(patientID));
            [rEmp, ~, nValid, eligible, reason] = EvaluatePair( ...
                entry, metabA, metabB, minValidParts, false);
            if prepareShifts
                shiftValues = ComputeCircularShiftPool(entry.partTable, metabA, metabB, minValidParts);
            else
                shiftValues = [];
            end
            shiftEligible = eligible && ~isempty(shiftValues) && any(isfinite(shiftValues));
            newRow = table(patientID, metabA, metabB, nValid, eligible, ...
                shiftEligible, rEmp, reason, ...
                'VariableNames', {'patientID','metaboliteA','metaboliteB', ...
                'nValidParts','eligible','shiftEligible','rEmpirical','reason'});
            newRow.shiftRValues = {shiftValues(:)};
            patientTable = [patientTable; newRow]; %#ok<AGROW>
        end
        rows = firstRow:height(patientTable);
        pairTable = [pairTable; table(metabA, metabB, ...
            sum(patientTable.eligible(rows)), sum(patientTable.shiftEligible(rows)), ...
            'VariableNames', {'metaboliteA','metaboliteB', ...
            'nEligiblePatients','nShiftEligiblePatients'})]; %#ok<AGROW>
    end
end

end

function values = ComputeCircularShiftPool(T, metabA, metabB, minValidParts)
values = [];
try
    if ismember("part", string(T.Properties.VariableNames))
        [~, order] = sort(T.part);
        T = T(order, :);
    end
    x = GetTableColumn(T, metabA);
    y = GetTableColumn(T, metabB);
catch
    return;
end
if sum(isfinite(x) & isfinite(y)) < minValidParts
    return;
end

values = nan(numel(x)-1, 1);
for shift = 1:numel(x)-1
    shifted = circshift(y, shift);
    valid = isfinite(x) & isfinite(shifted);
    if sum(valid) >= minValidParts
        values(shift) = corr(x(valid), shifted(valid));
    end
end
values = values(isfinite(values));
end

% -------------------------------------------------------------------------
function [rEmp, rModel, nValid, eligible, reason] = EvaluatePair(entry, metabA, metabB, minValidParts, requireModel)
rEmp = NaN;
rModel = NaN;
nValid = 0;
eligible = false;
reason = "";

try
    x = GetTableColumn(entry.partTable, metabA);
    y = GetTableColumn(entry.partTable, metabB);
catch
    reason = "missing empirical metabolite column";
    return;
end

valid = isfinite(x) & isfinite(y);
nValid = sum(valid);
if nValid < minValidParts
    reason = "insufficient valid parts";
    return;
end

rEmp = PearsonR(x(valid), y(valid));
if ~isfinite(rEmp)
    reason = "nonfinite empirical correlation";
    return;
end

if requireModel
    try
        rModel = GetSquareValue(entry.modelCorrelationTable, metabA, metabB);
    catch
        reason = "missing model correlation";
        return;
    end
    if ~isfinite(rModel)
        reason = "nonfinite model correlation";
        return;
    end
end

eligible = true;
end

% -------------------------------------------------------------------------
function T = BuildGlobalCRLBMajorityTable(modelOutputs, patientIDs, threshold)
allMetabs = strings(0, 1);
allValues = cell(0, 1);
for p = 1:numel(patientIDs)
    field = char(patientIDs(p));
    patient = modelOutputs.patientResultsByID.(field);
    if ~isfield(patient, 'partCRLB') || ~isfield(patient, 'metabList')
        error('Patient %s lacks partCRLB/metabList.', patientIDs(p));
    end
    metabs = string(patient.metabList(:));
    crlb = double(patient.partCRLB);
    for m = 1:numel(metabs)
        idx = FindMatchingNameIndex(allMetabs, metabs(m));
        vals = crlb(:, m);
        vals = vals(isfinite(vals));
        if isempty(idx)
            allMetabs(end+1, 1) = metabs(m); %#ok<AGROW>
            allValues{end+1, 1} = vals(:); %#ok<AGROW>
        else
            allValues{idx} = [allValues{idx}; vals(:)];
        end
    end
end

T = table();
for m = 1:numel(allMetabs)
    vals = allValues{m};
    nFinite = numel(vals);
    nUnder = sum(vals < threshold);
    nOver = sum(vals >= threshold);
    T = [T; table(allMetabs(m), nFinite, nUnder, nOver, ...
        100*nUnder/max(nFinite, 1), threshold, nUnder > nOver, ...
        'VariableNames', {'metabolite','nFiniteCRLB','nCRLBUnderThreshold', ...
        'nCRLBOverOrEqualThreshold','percentUnderThreshold','threshold','keepByCRLB'})]; %#ok<AGROW>
end
end

% -------------------------------------------------------------------------
function qualityTable = BuildTemporalCRLBQualityTable(dataByID, patientIDs, metabs, division, threshold, requiredFraction)
metabs = string(metabs(:));
goodCount = zeros(numel(metabs), 1);
instanceCount = zeros(numel(metabs), 1);

for p = 1:numel(patientIDs)
    entry = dataByID.(char(patientIDs(p)));
    coordTable = entry.empiricalCoordTable;
    partTable = entry.partTable;
    if isempty(coordTable)
        continue;
    end
    coordTable.name = string(coordTable.name);
    coordTable.filename = string(coordTable.filename);
    crlbColumn = FindCRLBColumn(coordTable);
    parsedDivision = nan(height(coordTable), 1);
    parsedPart = nan(height(coordTable), 1);
    for r = 1:height(coordTable)
        [~, baseName, ext] = fileparts(coordTable.filename(r));
        token = regexp(string(baseName) + string(ext), ...
            '.*Division_(\d+)_(?:part_)?(\d+)\.basis\.coord$', 'tokens', 'once');
        if ~isempty(token)
            parsedDivision(r) = str2double(token{1});
            parsedPart(r) = str2double(token{2});
        end
    end
    coordTable.filterDivision = parsedDivision;
    coordTable.filterPart = parsedPart;
    coordTable = coordTable(coordTable.filterDivision == division, :);
    parts = partTable.part;

    for m = 1:numel(metabs)
        values = nan(numel(parts), 1);
        for partIdx = 1:numel(parts)
            idx = coordTable.name == metabs(m) & coordTable.filterPart == parts(partIdx);
            if any(idx)
                tmp = coordTable.(char(crlbColumn))(idx);
                values(partIdx) = double(tmp(1));
            end
        end
        instanceCount(m) = instanceCount(m) + numel(values);
        goodCount(m) = goodCount(m) + sum(values < threshold, 'omitnan');
    end
end

fraction = goodCount ./ instanceCount;
fails = fraction < requiredFraction;
passes = ~fails & isfinite(fraction);
qualityTable = table(metabs, goodCount, instanceCount, fraction, fails, passes, ...
    'VariableNames', {'metabolite','nCRLBUnder100','nInstances', ...
    'fractionCRLBUnder100','fails90PercentRule','passesRequiredFraction'});
end

function name = FindCRLBColumn(T)
candidates = ["CRLB", "crlb", "SD", "sd", "percentSD", "PercentSD", ...
    "pctSD", "pctCrLB", "crlbPercent", "CRLBPercent"];
vars = string(T.Properties.VariableNames);
for k = 1:numel(candidates)
    if ismember(candidates(k), vars)
        name = candidates(k);
        return;
    end
end
error('Could not find a CRLB / %%SD column in empirical coord table.');
end

% -------------------------------------------------------------------------
function values = GetTableColumn(T, label)
idx = FindMatchingNameIndex(string(T.Properties.VariableNames), label);
if isempty(idx)
    error('Missing table column %s.', label);
end
values = double(T{:, idx});
values = values(:);
end

function value = GetSquareValue(T, rowLabel, colLabel)
if ~istable(T) || isempty(T)
    error('Model correlation table is unavailable.');
end
if ~isempty(T.Properties.RowNames)
    rows = string(T.Properties.RowNames);
else
    rows = string(T.Properties.VariableNames);
end
cols = string(T.Properties.VariableNames);
r = FindMatchingNameIndex(rows, rowLabel);
c = FindMatchingNameIndex(cols, colLabel);
if isempty(r) || isempty(c)
    error('Missing model correlation labels.');
end
value = double(T{r, c});
end

function tf = HasMatchingTableVariable(T, label)
tf = istable(T) && ~isempty(FindMatchingNameIndex(string(T.Properties.VariableNames), label));
end

function tf = HasMatchingSquareLabel(T, label)
if ~istable(T) || isempty(T)
    tf = false;
    return;
end
if ~isempty(T.Properties.RowNames)
    rows = string(T.Properties.RowNames);
else
    rows = string(T.Properties.VariableNames);
end
tf = ~isempty(FindMatchingNameIndex(rows, label)) && ...
    ~isempty(FindMatchingNameIndex(string(T.Properties.VariableNames), label));
end

function r = PearsonR(x, y)
x = double(x(:));
y = double(y(:));
if numel(x) < 3
    r = NaN;
    return;
end
x = x - mean(x);
y = y - mean(y);
denom = sqrt(sum(x.^2) * sum(y.^2));
if denom <= 0 || ~isfinite(denom)
    r = NaN;
else
    r = sum(x .* y) / denom;
end
end

function out = KeepMatchingNames(input, allowed)
input = string(input(:));
allowed = string(allowed(:));
keep = false(numel(input), 1);
for k = 1:numel(input)
    keep(k) = ~isempty(FindMatchingNameIndex(allowed, input(k)));
end
out = input(keep);
end

function idx = FindMatchingNameIndex(names, label)
names = string(names(:));
label = string(label);
idx = find(strcmp(names, label), 1, 'first');
if ~isempty(idx), return; end
labelValid = string(matlab.lang.makeValidName(char(label)));
validNames = strings(numel(names), 1);
for k = 1:numel(names)
    validNames(k) = string(matlab.lang.makeValidName(char(names(k))));
end
idx = find(strcmp(validNames, labelValid), 1, 'first');
end

function tf = NamesMatch(a, b)
tf = ~isempty(FindMatchingNameIndex(string(a), string(b)));
end

function value = GetOption(s, name, defaultValue)
if isstruct(s) && isfield(s, name)
    value = s.(name);
else
    value = defaultValue;
end
end

function value = GetOptionalField(s, name, defaultValue)
if isstruct(s) && isfield(s, name)
    value = s.(name);
else
    value = defaultValue;
end
end
