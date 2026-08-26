function outputs = CompareBasisCRLBToCurrentModels(covOutputs, deGraafOutputs, cfg, diagCfg)
%CompareBasisCRLBToCurrentModels Diagnostic basis-geometry model comparison.
%
% This function is deliberately separate from the production pairwise and
% Wishart tests. It reads the already-loaded LCModel fitted component traces,
% forms H = real(B' * B) in the active independent-component space, propagates
% the resulting covariance geometry to configured reported sums, and compares
% that one matrix with the two current patient-level model summaries.
%
% No baseline column is added. No production matrix is changed or replaced.

    if nargin < 3 || isempty(cfg)
        cfg = ProjectConfig();
    end
    if nargin < 4 || isempty(diagCfg)
        diagCfg = struct();
    end

    opts = ParseOptions(cfg, diagCfg);
    patientIDs = SelectPatientIDs(covOutputs, deGraafOutputs, opts.patientIDs);
    nPatients = numel(patientIDs);

    patientResults = struct([]);
    patientResultsByID = struct;

    fprintf('\nBasis-function covariance/CRLB geometry diagnostic\n');
    fprintf('Division: %d\n', opts.division);
    fprintf('Patients selected: %d\n', nPatients);
    fprintf('Baseline included: 0\n');
    fprintf(['Covariance scale: unresolved; normalized and bestFitScale ', ...
        'diagnostics only\n']);

    for pIdx = 1:nPatients
        patientID = patientIDs(pIdx);
        fprintf('Basis diagnostic %s (%d/%d)\n', patientID, pIdx, nPatients);

        covPatient = GetPatientResult(covOutputs, patientID);
        deGraafPatient = GetPatientResult(deGraafOutputs, patientID);

        patientResult = AnalyzePatient( ...
            patientID, covPatient, deGraafPatient, opts);

        if pIdx == 1
            patientResults = patientResult;
        else
            patientResults(pIdx) = patientResult;
        end
        safeID = matlab.lang.makeValidName(char(patientID));
        patientResultsByID.(safeID) = patientResult;
    end

    cohortTable = BuildCohortTable(patientResults);
    cohortSummaryTable = BuildCohortSummaryTable(cohortTable);
    group = BuildGroupComparison( ...
        covOutputs, deGraafOutputs, patientResults, patientIDs, opts);

    outputs = struct();
    outputs.settings = opts;
    outputs.patientIDs = patientIDs;
    outputs.patientResults = patientResults;
    outputs.patientResultsByID = patientResultsByID;
    outputs.cohortTable = cohortTable;
    outputs.cohortSummaryTable = cohortSummaryTable;
    outputs.basisInspection = BuildBasisInspection(patientResults, opts);
    outputs.group = group;
    outputs.scaleHandling = struct( ...
        'hasJustifiedSigmaSquared', false, ...
        'sigmaSquaredSource', "none", ...
        'description', ...
        "The .coord files contain a spectral residual, but the existing " + ...
        "estimator does not use a residual/noise scale and the exported " + ...
        "frequency-domain residual is not documented as whitened noise. " + ...
        "CGeometry is therefore not called an absolute CRLB covariance. " + ...
        "Covariance comparisons use trace normalization and a separately " + ...
        "labelled bestFitScale diagnostic.");

    if isfield(patientResultsByID, matlab.lang.makeValidName(char(opts.detailedPatientID)))
        detailField = matlab.lang.makeValidName(char(opts.detailedPatientID));
        outputs.(detailField) = patientResultsByID.(detailField);
        PrintDetailedPatientReport(patientResultsByID.(detailField));

        if opts.makeFigures
            [figureFiles, figureMetadata] = MakeDiagnosticFigures( ...
                patientResultsByID.(detailField), opts);
            outputs.(detailField).figureFiles = figureFiles;
            outputs.(detailField).figureMetadata = figureMetadata;
            outputs.patientResultsByID.(detailField).figureFiles = figureFiles;
            outputs.patientResultsByID.(detailField).figureMetadata = figureMetadata;

            idx = find(patientIDs == opts.detailedPatientID, 1, 'first');
            if ~isempty(idx)
                outputs.patientResults(idx).figureFiles = figureFiles;
                outputs.patientResults(idx).figureMetadata = figureMetadata;
            end
        end
    else
        warning('Detailed patient %s was not in the selected patient set.', ...
            opts.detailedPatientID);
    end

    if opts.makeGroupFigures
        [groupFigureFiles, groupFigureMetadata] = MakeGroupDiagnosticFigures( ...
            outputs.group, opts);
        outputs.group.figureFiles = groupFigureFiles;
        outputs.group.figureMetadata = groupFigureMetadata;
    end

    fprintf('\nCohort basis diagnostic summary\n');
    disp(cohortSummaryTable)
    PrintGroupComparisonReport(outputs.group);
end


function opts = ParseOptions(cfg, diagCfg)
    opts = struct();
    opts.patientIDs = "all";
    opts.detailedPatientID = "P01";
    opts.division = 1;
    opts.reportedNames = DefaultReportedNames();
    opts.finalPanel = DefaultFinalPanel();
    opts.independentCandidateNames = DefaultIndependentCandidateNames();
    opts.makeFigures = true;
    opts.makeGroupFigures = true;
    opts.saveFigures = true;
    opts.closeFigures = true;
    opts.figureVisible = "off";

    rootDir = string(pwd);
    if HasFieldOrProperty(cfg, "paths")
        paths = GetFieldOrProperty(cfg, "paths");
        if HasFieldOrProperty(paths, "rootDir")
            rootDir = string(GetFieldOrProperty(paths, "rootDir"));
        end
    end
    opts.outputDir = fullfile(rootDir, "BasisCRLBDiagnostics");

    optionNames = string(fieldnames(opts));
    providedNames = string(fieldnames(diagCfg));
    for k = 1:numel(providedNames)
        name = providedNames(k);
        if ~ismember(name, optionNames)
            error('Unknown basis diagnostic option: %s', name);
        end
        opts.(name) = diagCfg.(name);
    end

    if ischar(opts.patientIDs)
        opts.patientIDs = string(opts.patientIDs);
    else
        opts.patientIDs = string(opts.patientIDs(:));
    end
    opts.detailedPatientID = string(opts.detailedPatientID);
    opts.reportedNames = string(opts.reportedNames(:));
    opts.finalPanel = string(opts.finalPanel(:));
    opts.independentCandidateNames = string(opts.independentCandidateNames(:));
    opts.outputDir = string(opts.outputDir);
    opts.figureVisible = string(opts.figureVisible);

    if numel(unique(opts.reportedNames)) ~= numel(opts.reportedNames)
        error('reportedNames must be unique.');
    end
    if ~all(ismember(opts.finalPanel, opts.reportedNames))
        error('Every finalPanel name must occur in reportedNames.');
    end
end


function patientIDs = SelectPatientIDs(covOutputs, deGraafOutputs, requested)
    covIDs = GetPatientIDs(covOutputs);
    deGraafIDs = GetPatientIDs(deGraafOutputs);
    commonIDs = intersect(covIDs, deGraafIDs, 'stable');

    if isscalar(requested) && strcmpi(requested, "all")
        patientIDs = commonIDs;
    else
        requested = string(requested(:));
        missing = requested(~ismember(requested, commonIDs));
        if ~isempty(missing)
            error('Requested patients unavailable in both inputs: %s', ...
                strjoin(missing, ', '));
        end
        patientIDs = requested;
    end
end


function ids = GetPatientIDs(outputs)
    if ~isfield(outputs, 'patientResults') || isempty(outputs.patientResults)
        ids = strings(0, 1);
        return;
    end
    ids = string({outputs.patientResults.patientID}).';
end


function patient = GetPatientResult(outputs, patientID)
    ids = GetPatientIDs(outputs);
    idx = find(ids == patientID, 1, 'first');
    if isempty(idx)
        error('Patient %s is not present in the supplied outputs.', patientID);
    end
    patient = outputs.patientResults(idx);
end


function result = AnalyzePatient(patientID, covPatient, deGraafPatient, opts)
    reportedNames = opts.reportedNames;
    nReported = numel(reportedNames);
    storeFullMatrices = patientID == opts.detailedPatientID;

    if double(deGraafPatient.division) ~= opts.division || ...
            double(covPatient.division) ~= opts.division
        error('Patient %s inputs do not both represent Division-%d.', ...
            patientID, opts.division);
    end

    partInfo = deGraafPatient.partInfo;
    usePart = partInfo.wasParsed & isfinite(partInfo.part);
    partInfo = partInfo(usePart, :);
    partInfo = partInfo(partInfo.part >= 1, :);
    [~, order] = sort(partInfo.part);
    partInfo = partInfo(order, :);

    nParts = height(partInfo);
    if nParts == 0
        error('No parsed Division-%d parts are available for %s.', ...
            opts.division, patientID);
    end

    CReportedStack = nan(nReported, nReported, nParts);
    partDiagnostics = repmat(EmptyPartDiagnostic(), nParts, 1);
    componentNamesByPart = cell(nParts, 1);
    excludedNamesByPart = cell(nParts, 1);
    allActualNamesByPart = cell(nParts, 1);

    for partIdx = 1:nParts
        coordFile = string(partInfo.coordFile(partIdx));
        [fitData, quantRows] = GetLoadedPartData(covPatient, coordFile);

        [B, componentNames, construction] = BuildIndependentBasis( ...
            fitData, quantRows);
        allActualNamesByPart{partIdx} = construction.actualComponentNames;
        componentNamesByPart{partIdx} = componentNames;

        excluded = setdiff(opts.independentCandidateNames, componentNames, 'stable');
        excludedNamesByPart{partIdx} = excluded;

        [CComponent, hDiag] = SolveGeometry(B);
        [M, propagationTable] = BuildReportedTransformation( ...
            componentNames, reportedNames, opts.independentCandidateNames);

        if hDiag.numericallyUsable
            CReported = M * CComponent * M.';
            CReported = (CReported + CReported.') ./ 2;
            CReportedStack(:, :, partIdx) = CReported;
        end

        partDiag = EmptyPartDiagnostic();
        partDiag.part = partInfo.part(partIdx);
        partDiag.coordFile = coordFile;
        partDiag.spectralPointCount = size(B, 1);
        partDiag.axisStartPPM = FirstOrNaN(fitData.axis);
        partDiag.axisEndPPM = LastOrNaN(fitData.axis);
        partDiag.isReal = isreal(B);
        partDiag.actualComponentNames = construction.actualComponentNames;
        partDiag.componentNames = componentNames;
        partDiag.excludedCandidateNames = excluded;
        partDiag.exclusionTable = construction.exclusionTable;
        partDiag.componentConcentrations = construction.componentConcentrations;
        if storeFullMatrices
            partDiag.B = B;
            partDiag.H = hDiag.H;
            partDiag.CComponentGeometry = CComponent;
            partDiag.M = M;
        end
        partDiag.propagationTable = propagationTable;
        partDiag.CBasisReported = CReportedStack(:, :, partIdx);
        partDiag.sizeH = hDiag.sizeH;
        partDiag.rankH = hDiag.rankH;
        partDiag.conditionNumber = hDiag.conditionNumber;
        partDiag.minimumEigenvalue = hDiag.minimumEigenvalue;
        partDiag.maximumEigenvalue = hDiag.maximumEigenvalue;
        partDiag.positiveDefinite = hDiag.positiveDefinite;
        partDiag.numericallyUsable = hDiag.numericallyUsable;
        partDiag.rcondH = hDiag.rcondH;
        partDiag.residualVarianceObservedNotUsed = ...
            ObservedResidualVariance(fitData);
        partDiagnostics(partIdx) = partDiag;
    end

    usableParts = [partDiagnostics.numericallyUsable].';
    if ~any(usableParts)
        error('Every basis information matrix was unusable for %s.', patientID);
    end

    CBasisReported = mean(CReportedStack(:, :, usableParts), 3, 'omitnan');
    CBasisReported = (CBasisReported + CBasisReported.') ./ 2;
    RBasisReported = CovarianceToCorrelation(CBasisReported);

    RDeGraaf = ReorderMetaboliteTable( ...
        deGraafPatient.meanAmplitudeCorrTable, reportedNames);
    CWishartCurrent = ReorderMetaboliteTable( ...
        deGraafPatient.meanAmplitudeCovTable, reportedNames);
    RFromCurrentWishart = CovarianceToCorrelation(CWishartCurrent);

    validCorrelationNames = ValidCorrelationNames( ...
        reportedNames, RBasisReported, RDeGraaf);
    validCovarianceNames = ValidCovarianceNames( ...
        reportedNames, CBasisReported, CWishartCurrent);
    validCurrentMethodNames = ValidCorrelationNames( ...
        reportedNames, RFromCurrentWishart, RDeGraaf);

    corrConfigured = CompareCorrelation( ...
        RBasisReported, RDeGraaf, reportedNames, validCorrelationNames, ...
        "R_basis", "R_DeGraaf");
    corrFinal = CompareCorrelation( ...
        RBasisReported, RDeGraaf, reportedNames, ...
        intersect(opts.finalPanel, validCorrelationNames, 'stable'), ...
        "R_basis", "R_DeGraaf");

    covarianceConfigured = CompareCovariance( ...
        CBasisReported, CWishartCurrent, reportedNames, validCovarianceNames);
    covarianceFinal = CompareCovariance( ...
        CBasisReported, CWishartCurrent, reportedNames, ...
        intersect(opts.finalPanel, validCovarianceNames, 'stable'));

    currentConfigured = CompareCorrelation( ...
        RFromCurrentWishart, RDeGraaf, reportedNames, ...
        validCurrentMethodNames, "R_from_CurrentWishart", "R_DeGraaf");
    currentFinal = CompareCorrelation( ...
        RFromCurrentWishart, RDeGraaf, reportedNames, ...
        intersect(opts.finalPanel, validCurrentMethodNames, 'stable'), ...
        "R_from_CurrentWishart", "R_DeGraaf");

    componentFrequencyTable = BuildComponentFrequencyTable( ...
        opts.independentCandidateNames, componentNamesByPart, ...
        allActualNamesByPart, nParts);
    hSummaryTable = BuildHSummaryTable(partDiagnostics);

    result = struct();
    result.patientID = patientID;
    result.componentNames = unique(vertcat(componentNamesByPart{:}), 'stable');
    result.reportedNames = reportedNames;
    result.independentCandidateNames = opts.independentCandidateNames;
    result.componentNamesByPart = componentNamesByPart;
    result.excludedNamesByPart = excludedNamesByPart;
    result.componentFrequencyTable = componentFrequencyTable;
    result.partDiagnostics = partDiagnostics;
    result.H = {partDiagnostics.H}.';
    result.CComponentGeometry = {partDiagnostics.CComponentGeometry}.';
    result.CBasisReportedByPart = CReportedStack;
    result.CBasisReported = CBasisReported;
    result.RBasisReported = RBasisReported;
    result.CWishartCurrent = CWishartCurrent;
    result.RDeGraaf = RDeGraaf;
    result.RFromCurrentWishart = RFromCurrentWishart;
    result.validCorrelationNames = validCorrelationNames;
    result.validCovarianceNames = validCovarianceNames;
    result.validCurrentMethodNames = validCurrentMethodNames;
    result.comparisons.basisVsDeGraaf.configured = corrConfigured;
    result.comparisons.basisVsDeGraaf.final11 = corrFinal;
    result.comparisons.basisVsWishart.configured = covarianceConfigured;
    result.comparisons.basisVsWishart.final11 = covarianceFinal;
    result.comparisons.currentWishartRVsDeGraaf.configured = currentConfigured;
    result.comparisons.currentWishartRVsDeGraaf.final11 = currentFinal;
    result.hSummaryTable = hSummaryTable;
    result.conditionNumberSummary = SummarizeVector(hSummaryTable.conditionNumber);
    result.nParts = nParts;
    result.nUsableParts = sum(usableParts);
    result.basisSource = ...
        "covOutputs.patientResults.fitData.basisData, parsed by " + ...
        "VDIIO.ReadLCMCoord from LCModel .coord individual fitted-component " + ...
        "plots; parser-subtracted LCModel baseline; no baseline column added";
    result.basisScaling = ...
        "Each exported fitted component trace was divided by its matching " + ...
        ".coord fitted concentration before H was formed. The exported trace " + ...
        "is concentration-scaled; division expresses B as the local spectral " + ...
        "derivative with respect to the reported concentration parameter.";
    result.noiseScale = ...
        "No justified sigma^2 found or used. Residual variance is retained " + ...
        "only as an observed diagnostic and is not treated as whitened noise.";
end


function d = EmptyPartDiagnostic()
    d = struct( ...
        'part', NaN, ...
        'coordFile', "", ...
        'spectralPointCount', NaN, ...
        'axisStartPPM', NaN, ...
        'axisEndPPM', NaN, ...
        'isReal', false, ...
        'actualComponentNames', strings(0, 1), ...
        'componentNames', strings(0, 1), ...
        'excludedCandidateNames', strings(0, 1), ...
        'exclusionTable', table(), ...
        'componentConcentrations', zeros(0, 1), ...
        'B', [], ...
        'H', [], ...
        'CComponentGeometry', [], ...
        'M', [], ...
        'propagationTable', table(), ...
        'CBasisReported', [], ...
        'sizeH', [0 0], ...
        'rankH', NaN, ...
        'conditionNumber', NaN, ...
        'minimumEigenvalue', NaN, ...
        'maximumEigenvalue', NaN, ...
        'positiveDefinite', false, ...
        'numericallyUsable', false, ...
        'rcondH', NaN, ...
        'residualVarianceObservedNotUsed', NaN);
end


function [fitData, quantRows] = GetLoadedPartData(covPatient, coordFile)
    coordFiles = string(covPatient.coordFiles(:));
    targetFull = NormalizePath(coordFile);
    normalized = lower(replace(coordFiles, '/', '\'));
    idx = find(normalized == targetFull, 1, 'first');

    if isempty(idx)
        targetName = string(GetFileName(coordFile));
        names = arrayfun(@GetFileName, coordFiles);
        idx = find(strcmpi(names, targetName), 1, 'first');
    end
    if isempty(idx)
        error('Loaded fitData does not contain %s.', coordFile);
    end

    fitData = covPatient.fitData(idx);
    T = covPatient.coordTable;
    targetName = string(GetFileName(coordFile));
    rowMask = strcmpi(string(T.filename), targetName);
    quantRows = T(rowMask, :);
    if isempty(quantRows)
        error('Loaded quantification table has no rows for %s.', targetName);
    end
end


function p = NormalizePath(p)
    p = replace(string(p), '/', '\');
    p = lower(p);
end


function name = GetFileName(pathValue)
    [~, base, ext] = fileparts(char(pathValue));
    name = string([base, ext]);
end


function [B, componentNames, info] = BuildIndependentBasis(fitData, quantRows)
    rawB = double(fitData.basisData);
    rawNames = string(fitData.basisMetName(:));

    if isempty(rawB) || size(rawB, 2) ~= numel(rawNames)
        error('fitData basisData/basisMetName dimensions are inconsistent.');
    end

    actualNames = strings(numel(rawNames), 1);
    concentration = nan(numel(rawNames), 1);
    keep = true(numel(rawNames), 1);
    reasons = strings(numel(rawNames), 1);

    quantNames = string(quantRows.name);
    for k = 1:numel(rawNames)
        [actualNames(k), qIdx] = ResolveComponentName(rawNames(k), quantNames);

        if contains(actualNames(k), "+")
            keep(k) = false;
            reasons(k) = "reported sum is derived, not independent";
            continue;
        end
        if isempty(qIdx)
            keep(k) = false;
            reasons(k) = "no matching fitted concentration";
            continue;
        end

        concentration(k) = double(quantRows.sig(qIdx));
        if ~isfinite(concentration(k)) || concentration(k) == 0
            keep(k) = false;
            reasons(k) = "zero or nonfinite fitted concentration";
        elseif any(~isfinite(rawB(:, k)))
            keep(k) = false;
            reasons(k) = "nonfinite spectral samples";
        elseif norm(rawB(:, k)) == 0
            keep(k) = false;
            reasons(k) = "zero spectral column";
        end
    end

    if numel(unique(actualNames(keep))) ~= sum(keep)
        error('Duplicate active component names occur in fitted basis data.');
    end

    componentNames = actualNames(keep);
    componentConcentrations = concentration(keep);
    B = rawB(:, keep) ./ componentConcentrations.';

    excludedMask = ~keep;
    info = struct();
    info.actualComponentNames = actualNames;
    info.componentConcentrations = componentConcentrations;
    info.exclusionTable = table( ...
        rawNames(excludedMask), actualNames(excludedMask), reasons(excludedMask), ...
        'VariableNames', {'parsedName', 'resolvedName', 'reason'});
end


function [resolved, idx] = ResolveComponentName(parsedName, quantNames)
    parsedName = string(parsedName);
    idx = find(strcmpi(quantNames, parsedName), 1, 'first');
    resolved = parsedName;

    % VDIIO's current basis-header regular expression intentionally captures
    % word characters, so LCModel's leading minus in -CrCH2 is not retained.
    if isempty(idx) && strcmpi(parsedName, "CrCH2")
        idx = find(strcmpi(quantNames, "-CrCH2"), 1, 'first');
        if ~isempty(idx)
            resolved = "-CrCH2";
        end
    end
end


function [C, d] = SolveGeometry(B)
    H = real(B' * B);
    H = (H + H.') ./ 2;
    n = size(H, 1);

    d = struct();
    d.H = H;
    d.sizeH = size(H);
    d.rankH = rank(H);
    d.conditionNumber = cond(H);
    d.rcondH = rcond(H);

    eigenvalues = eig(H, 'vector');
    d.minimumEigenvalue = min(real(eigenvalues));
    d.maximumEigenvalue = max(real(eigenvalues));
    [~, cholStatus] = chol(H);
    d.positiveDefinite = cholStatus == 0;

    tolerance = max(1, n) * eps(max(1, norm(H, 1)));
    d.numericallyUsable = d.positiveDefinite && ...
        d.rankH == n && isfinite(d.conditionNumber) && ...
        d.rcondH > tolerance / max(1, norm(H, 1));

    if ~d.numericallyUsable
        C = nan(n);
        warning(['Basis H is singular or numerically unusable: size %dx%d, ', ...
            'rank %d, cond %.6g, minEig %.6g. No pseudoinverse was used.'], ...
            n, n, d.rankH, d.conditionNumber, d.minimumEigenvalue);
        return;
    end

    C = H \ eye(n);
    C = (C + C.') ./ 2;
end


function [M, propagationTable] = BuildReportedTransformation( ...
        componentNames, reportedNames, candidateNames)
    nReported = numel(reportedNames);
    nComponents = numel(componentNames);
    M = zeros(nReported, nComponents);
    termsText = strings(nReported, 1);
    activeTermsText = strings(nReported, 1);
    inactiveTermsText = strings(nReported, 1);

    for r = 1:nReported
        terms = split(reportedNames(r), "+");
        terms = string(terms(:));
        termsText(r) = strjoin(terms, " + ");
        activeTerms = strings(0, 1);
        inactiveTerms = strings(0, 1);

        for t = 1:numel(terms)
            idx = find(componentNames == terms(t), 1, 'first');
            if ~isempty(idx)
                M(r, idx) = M(r, idx) + 1;
                activeTerms(end+1, 1) = terms(t); %#ok<AGROW>
            else
                if ~ismember(terms(t), candidateNames)
                    error('Reported term %s is not a known independent candidate.', terms(t));
                end
                inactiveTerms(end+1, 1) = terms(t); %#ok<AGROW>
            end
        end
        activeTermsText(r) = strjoin(activeTerms, " + ");
        inactiveTermsText(r) = strjoin(inactiveTerms, " + ");
    end

    propagationTable = table( ...
        reportedNames, termsText, activeTermsText, inactiveTermsText, ...
        sum(M ~= 0, 2), ...
        'VariableNames', {'reportedName', 'definedTerms', 'activeTerms', ...
        'inactiveFixedTerms', 'nActiveTerms'});
end


function value = ObservedResidualVariance(fitData)
    residual = real(double(fitData.residual(:)));
    residual = residual(isfinite(residual));
    if numel(residual) < 2
        value = NaN;
    else
        value = var(residual, 1);
    end
end


function value = FirstOrNaN(x)
    if isempty(x), value = NaN; else, value = double(x(1)); end
end


function value = LastOrNaN(x)
    if isempty(x), value = NaN; else, value = double(x(end)); end
end


function A = ReorderMetaboliteTable(T, names)
    rowNames = string(T.Properties.RowNames);
    columnNames = string(T.Properties.VariableNames);
    requestedColumnNames = string(matlab.lang.makeValidName(cellstr(names)));
    [rowFound, rowIdx] = ismember(names, rowNames);
    [columnFound, columnIdx] = ismember(requestedColumnNames, columnNames);

    A = nan(numel(names));
    found = rowFound & columnFound;
    A(found, found) = table2array(T(rowIdx(found), columnIdx(found)));
    A = double(A);
end


function R = CovarianceToCorrelation(C)
    C = double(C);
    C = (C + C.') ./ 2;
    variance = diag(C);
    denom = sqrt(variance * variance.');
    R = C ./ denom;
    invalid = ~isfinite(variance) | variance <= 0;
    R(invalid, :) = NaN;
    R(:, invalid) = NaN;
    valid = ~invalid;
    diagonalIndices = find(valid);
    R(sub2ind(size(R), diagonalIndices, diagonalIndices)) = 1;
    R = (R + R.') ./ 2;

    finiteMask = isfinite(R);
    tolerance = 1e-10;
    if any(abs(R(finiteMask)) > 1 + tolerance)
        warning('A covariance-derived correlation exceeded [-1,1] beyond tolerance.');
    end
    R(finiteMask) = min(1, max(-1, R(finiteMask)));
end


function names = ValidCorrelationNames(allNames, A, B)
    valid = isfinite(diag(A)) & isfinite(diag(B)) & ...
        diag(A) > 0 & diag(B) > 0;
    names = allNames(valid);
end


function names = ValidCovarianceNames(allNames, A, B)
    valid = isfinite(diag(A)) & isfinite(diag(B)) & ...
        diag(A) > 0 & diag(B) > 0;
    names = allNames(valid);
end


function comparison = CompareCorrelation(A, B, allNames, selectedNames, aLabel, bLabel)
    idx = FindNameIndices(allNames, selectedNames);
    A = A(idx, idx);
    B = B(idx, idx);
    [x, y, nameA, nameB] = UniqueOffDiagonalEntries(B, A, selectedNames);
    finite = isfinite(x) & isfinite(y);
    x = x(finite);
    y = y(finite);
    nameA = nameA(finite);
    nameB = nameB(finite);

    difference = y - x;
    comparison = struct();
    comparison.names = selectedNames;
    comparison.nMetabolites = numel(selectedNames);
    comparison.nPairs = numel(x);
    comparison.matrixCorrelation = Pearson(x, y);
    comparison.RMSE = RootMeanSquare(difference);
    comparison.MAD = MeanAbsolute(difference);
    comparison.medianAbsoluteDifference = MedianAbsolute(difference);
    comparison.meanSignedDifference = MeanFinite(difference);
    comparison.maximumAbsoluteDifference = MaxAbsolute(difference);
    comparison.signAgreement = MeanFinite(double(sign(x) == sign(y)));
    comparison.referenceLabel = string(bLabel);
    comparison.comparisonLabel = string(aLabel);
    comparison.referenceEntries = x;
    comparison.comparisonEntries = y;
    comparison.topDiscrepancies = BuildTopDiscrepancyTable( ...
        nameA, nameB, y, x, aLabel, bLabel);
end


function comparison = CompareCovariance(CBasis, CWishart, allNames, selectedNames)
    idx = FindNameIndices(allNames, selectedNames);
    A = CBasis(idx, idx);
    B = CWishart(idx, idx);

    traceA = trace(A);
    traceB = trace(B);
    comparison = struct();
    comparison.names = selectedNames;
    comparison.nMetabolites = numel(selectedNames);
    comparison.normalization = "C / trace(C)";
    comparison.hasJustifiedSigmaSquared = false;

    if ~isfinite(traceA) || ~isfinite(traceB) || traceA <= 0 || traceB <= 0
        comparison.normalizedBasis = nan(size(A));
        comparison.normalizedWishart = nan(size(B));
        comparison.nEntries = 0;
        comparison.matrixCorrelation = NaN;
        comparison.RMSE = NaN;
        comparison.MAD = NaN;
        comparison.diagonalVarianceCorrelation = NaN;
        comparison.bestFitScale = NaN;
        comparison.bestScaledRMSE = NaN;
        comparison.bestScaledMAD = NaN;
        return;
    end

    Anorm = A ./ traceA;
    Bnorm = B ./ traceB;
    mask = triu(true(size(A)));
    finite = mask & isfinite(Anorm) & isfinite(Bnorm);
    x = Bnorm(finite);
    y = Anorm(finite);
    difference = y - x;

    comparison.normalizedBasis = Anorm;
    comparison.normalizedWishart = Bnorm;
    comparison.nEntries = numel(x);
    comparison.matrixCorrelation = Pearson(x, y);
    comparison.RMSE = RootMeanSquare(difference);
    comparison.MAD = MeanAbsolute(difference);
    comparison.diagonalVarianceCorrelation = Pearson(diag(Bnorm), diag(Anorm));

    fullFinite = isfinite(A) & isfinite(B);
    denom = sum(A(fullFinite).^2);
    if denom > 0
        bestFitScale = sum(A(fullFinite) .* B(fullFinite)) ./ denom;
    else
        bestFitScale = NaN;
    end
    comparison.bestFitScale = bestFitScale;

    bestScaledDifference = bestFitScale .* A(mask) - B(mask);
    bestScaledDifference = bestScaledDifference(isfinite(bestScaledDifference));
    comparison.bestScaledRMSE = RootMeanSquare(bestScaledDifference);
    comparison.bestScaledMAD = MeanAbsolute(bestScaledDifference);
end


function idx = FindNameIndices(allNames, selectedNames)
    [found, idx] = ismember(selectedNames, allNames);
    if ~all(found)
        error('Selected names are not all present in the configured ordering.');
    end
end


function [x, y, nameA, nameB] = UniqueOffDiagonalEntries(A, B, names)
    n = numel(names);
    mask = triu(true(n), 1);
    [row, col] = find(mask);
    linear = sub2ind([n n], row, col);
    x = A(linear);
    y = B(linear);
    nameA = names(row);
    nameB = names(col);
end


function value = Pearson(x, y)
    finite = isfinite(x) & isfinite(y);
    x = double(x(finite));
    y = double(y(finite));
    if numel(x) < 2 || std(x) == 0 || std(y) == 0
        value = NaN;
        return;
    end
    R = corrcoef(x, y);
    value = R(1, 2);
end


function value = RootMeanSquare(x)
    x = x(isfinite(x));
    if isempty(x), value = NaN; else, value = sqrt(mean(x.^2)); end
end


function value = MeanAbsolute(x)
    x = x(isfinite(x));
    if isempty(x), value = NaN; else, value = mean(abs(x)); end
end


function value = MedianAbsolute(x)
    x = x(isfinite(x));
    if isempty(x), value = NaN; else, value = median(abs(x)); end
end


function value = MeanFinite(x)
    x = x(isfinite(x));
    if isempty(x), value = NaN; else, value = mean(x); end
end


function value = MaxAbsolute(x)
    x = x(isfinite(x));
    if isempty(x), value = NaN; else, value = max(abs(x)); end
end


function T = BuildTopDiscrepancyTable(nameA, nameB, comparisonValues, ...
        referenceValues, comparisonLabel, referenceLabel)
    difference = comparisonValues - referenceValues;
    absDifference = abs(difference);
    [~, order] = sort(absDifference, 'descend');
    order = order(1:min(10, numel(order)));

    comparisonVar = char(matlab.lang.makeValidName(comparisonLabel));
    referenceVar = char(matlab.lang.makeValidName(referenceLabel));
    T = table(nameA(order), nameB(order), comparisonValues(order), ...
        referenceValues(order), difference(order), absDifference(order), ...
        'VariableNames', {'metaboliteA', 'metaboliteB', comparisonVar, ...
        referenceVar, 'difference', 'absDifference'});
end


function T = BuildComponentFrequencyTable(candidateNames, usedByPart, ...
        actualByPart, nParts)
    actualUnion = unique(vertcat(actualByPart{:}), 'stable');
    names = unique([candidateNames; actualUnion], 'stable');
    nActual = zeros(numel(names), 1);
    nUsed = zeros(numel(names), 1);
    for k = 1:numel(names)
        nActual(k) = sum(cellfun(@(x) any(x == names(k)), actualByPart));
        nUsed(k) = sum(cellfun(@(x) any(x == names(k)), usedByPart));
    end
    T = table(names, nActual, nUsed, nParts - nUsed, ...
        ismember(names, candidateNames), ...
        'VariableNames', {'component', 'nPartsExported', 'nPartsUsed', ...
        'nPartsInactiveOrExcluded', 'inAudited28CandidateSet'});
end


function T = BuildHSummaryTable(partDiagnostics)
    part = reshape([partDiagnostics.part], [], 1);
    componentCount = reshape(arrayfun(@(x)x.sizeH(1), partDiagnostics), [], 1);
    rankH = reshape([partDiagnostics.rankH], [], 1);
    conditionNumber = reshape([partDiagnostics.conditionNumber], [], 1);
    rcondH = reshape([partDiagnostics.rcondH], [], 1);
    minimumEigenvalue = reshape([partDiagnostics.minimumEigenvalue], [], 1);
    maximumEigenvalue = reshape([partDiagnostics.maximumEigenvalue], [], 1);
    positiveDefinite = reshape([partDiagnostics.positiveDefinite], [], 1);
    numericallyUsable = reshape([partDiagnostics.numericallyUsable], [], 1);
    residualVarianceObservedNotUsed = reshape( ...
        [partDiagnostics.residualVarianceObservedNotUsed], [], 1);
    T = table( ...
        part, componentCount, rankH, conditionNumber, rcondH, ...
        minimumEigenvalue, maximumEigenvalue, positiveDefinite, ...
        numericallyUsable, residualVarianceObservedNotUsed, ...
        'VariableNames', {'part', 'componentCount', 'rankH', ...
        'conditionNumber', 'rcondH', 'minimumEigenvalue', ...
        'maximumEigenvalue', 'positiveDefinite', 'numericallyUsable', ...
        'residualVarianceObservedNotUsed'});
end


function s = SummarizeVector(x)
    x = double(x(:));
    x = x(isfinite(x));
    if isempty(x)
        s = struct('median', NaN, 'minimum', NaN, 'maximum', NaN, ...
            'Q1', NaN, 'Q3', NaN);
        return;
    end
    s = struct('median', median(x), 'minimum', min(x), 'maximum', max(x), ...
        'Q1', SimpleQuantile(x, 0.25), 'Q3', SimpleQuantile(x, 0.75));
end


function inspection = BuildBasisInspection(patientResults, opts)
    spectralPointCounts = [];
    componentCounts = [];
    allReal = true;
    activeSetVariesWithinPatients = false;
    actualUnion = strings(0, 1);

    for pIdx = 1:numel(patientResults)
        p = patientResults(pIdx);
        spectralPointCounts = [spectralPointCounts; ...
            reshape([p.partDiagnostics.spectralPointCount], [], 1)]; %#ok<AGROW>
        componentCounts = [componentCounts; ...
            p.hSummaryTable.componentCount]; %#ok<AGROW>
        allReal = allReal && all([p.partDiagnostics.isReal]);
        actualUnion = unique([actualUnion; p.componentNames], 'stable');

        referenceSet = sort(p.componentNamesByPart{1});
        for partIdx = 2:numel(p.componentNamesByPart)
            if ~isequal(sort(p.componentNamesByPart{partIdx}), referenceSet)
                activeSetVariesWithinPatients = true;
                break;
            end
        end
    end

    inspection = struct();
    inspection.source = [ ...
        "LCModel .coord fitted-component spectra already parsed into " + ...
        "covOutputs.patientResults.fitData.basisData"];
    inspection.domain = "real, phased frequency-domain spectra on the LCModel PPM axis";
    inspection.spectralPointCounts = unique(spectralPointCounts);
    inspection.componentCountRange = [min(componentCounts), max(componentCounts)];
    inspection.allMatricesReal = allReal;
    inspection.actualIndependentComponentUnion = actualUnion;
    inspection.audited28CandidateNames = opts.independentCandidateNames;
    inspection.unexpectedIndependentNames = setdiff( ...
        actualUnion, opts.independentCandidateNames, 'stable');
    inspection.candidateNamesNeverExported = setdiff( ...
        opts.independentCandidateNames, actualUnion, 'stable');
    inspection.activeSetVariesBetweenDivision1Parts = ...
        activeSetVariesWithinPatients;
    inspection.isPatientSpecific = numel(patientResults) > 1;
    inspection.patientSpecificExplanation = [ ...
        "The exported traces are LCModel fitted contributions after each " + ...
        "part's fitted shift, linewidth, phase, and amplitude. Active column " + ...
        "sets vary by part; the diagnostic therefore constructs B separately " + ...
        "for every patient and Division-1 part."];
end


function q = SimpleQuantile(x, probability)
    x = sort(double(x(:)));
    if isempty(x)
        q = NaN;
    elseif isscalar(x)
        q = x;
    else
        position = 1 + (numel(x) - 1) * probability;
        lowerIndex = floor(position);
        upperIndex = ceil(position);
        weight = position - lowerIndex;
        q = (1 - weight) * x(lowerIndex) + weight * x(upperIndex);
    end
end


function group = BuildGroupComparison(covOutputs, deGraafOutputs, ...
        patientResults, patientIDs, opts)
    empiricalTable = covOutputs.group.meanCorrTable;
    deGraafTable = deGraafOutputs.group.meanAmplitudeCorrTable;
    panelNames = string(empiricalTable.Properties.RowNames);
    deGraafNames = string(deGraafTable.Properties.RowNames);

    if ~isequal(panelNames, deGraafNames)
        error(['Existing empirical and DeGraaf group matrices do not have ', ...
            'the same metabolite ordering.']);
    end
    if ~all(ismember(panelNames, opts.reportedNames))
        error('The old group matrix panel is not available in the reported basis panel.');
    end

    % These two matrices are reused exactly from the objects plotted by
    % RunExperiment; neither is reconstructed from patient data here.
    REmpirical = double(table2array(empiricalTable));
    RDeGraaf = double(table2array(deGraafTable));

    nMetabolites = numel(panelNames);
    nPatients = numel(patientResults);
    basisPatientStack = nan(nMetabolites, nMetabolites, nPatients);
    for pIdx = 1:nPatients
        patient = patientResults(pIdx);
        idx = FindNameIndices(patient.reportedNames, panelNames);
        basisPatientStack(:, :, pIdx) = patient.RBasisReported(idx, idx);
    end

    [RBasis, basisPatientCountMatrix] = ...
        FisherAverageCorrelationStack(basisPatientStack);

    upperMask = triu(true(nMetabolites), 1);
    commonPairMask = upperMask & isfinite(REmpirical) & ...
        isfinite(RDeGraaf) & isfinite(RBasis);
    commonSymmetricMask = commonPairMask | commonPairMask.' | eye(nMetabolites) > 0;

    empiricalForComparison = REmpirical;
    deGraafForComparison = RDeGraaf;
    basisForComparison = RBasis;
    empiricalForComparison(~commonSymmetricMask) = NaN;
    deGraafForComparison(~commonSymmetricMask) = NaN;
    basisForComparison(~commonSymmetricMask) = NaN;

    empiricalVsDeGraaf = CompareCorrelation( ...
        empiricalForComparison, deGraafForComparison, panelNames, panelNames, ...
        "R_Empirical", "R_DeGraaf");
    empiricalVsBasis = CompareCorrelation( ...
        empiricalForComparison, basisForComparison, panelNames, panelNames, ...
        "R_Empirical", "R_Basis");
    basisVsDeGraaf = CompareCorrelation( ...
        basisForComparison, deGraafForComparison, panelNames, panelNames, ...
        "R_Basis", "R_DeGraaf");

    empiricalCountMatrix = ReorderCountMatrix( ...
        covOutputs.group.nPatientsCorrMatrix, covOutputs.metabList, panelNames);
    deGraafCountMatrix = ReorderCountMatrix( ...
        deGraafOutputs.group.nPatientsCorrMatrix, deGraafOutputs.metabList, panelNames);

    group = struct();
    group.metaboliteNames = panelNames;
    group.nMetabolites = nMetabolites;
    group.patientIDs = patientIDs;
    group.totalPatientsAvailable = nPatients;
    group.REmpirical = REmpirical;
    group.RDeGraaf = RDeGraaf;
    group.RBasis = RBasis;
    group.REmpiricalTable = empiricalTable;
    group.RDeGraafTable = deGraafTable;
    group.RBasisTable = MatrixToNamedTable(RBasis, panelNames);
    group.basisPatientStack = basisPatientStack;
    group.patientCountMatrix = basisPatientCountMatrix;
    group.patientCountMatrices.empirical = empiricalCountMatrix;
    group.patientCountMatrices.deGraaf = deGraafCountMatrix;
    group.patientCountMatrices.basis = basisPatientCountMatrix;
    group.commonPairMask = commonPairMask;
    group.nCommonOffDiagonalPairs = sum(commonPairMask, 'all');
    group.aggregation = [ ...
        "Elementwise Fisher average across patient RBasisReported matrices: " + ...
        "atanh after clamping to +/-0.999999, mean with missing patients " + ...
        "omitted, then tanh. This matches the existing empirical and " + ...
        "DeGraaf group correlation convention."];
    group.comparisons.empiricalVsDeGraaf = empiricalVsDeGraaf;
    group.comparisons.empiricalVsBasis = empiricalVsBasis;
    group.comparisons.basisVsDeGraaf = basisVsDeGraaf;
end


function [meanCorrelation, countMatrix] = FisherAverageCorrelationStack(stack)
    [nRows, nCols, ~] = size(stack);
    meanCorrelation = nan(nRows, nCols);
    countMatrix = sum(isfinite(stack), 3);

    for row = 1:nRows
        for column = 1:nCols
            values = squeeze(stack(row, column, :));
            values = values(isfinite(values));
            if isempty(values)
                continue;
            end
            if row == column
                meanCorrelation(row, column) = 1;
                continue;
            end
            values(values >= 1) = 0.999999;
            values(values <= -1) = -0.999999;
            meanCorrelation(row, column) = tanh(mean(atanh(values), 'omitnan'));
        end
    end
    meanCorrelation = (meanCorrelation + meanCorrelation.') ./ 2;
end


function matrixOut = ReorderCountMatrix(matrixIn, inputNames, requestedNames)
    inputNames = string(inputNames(:));
    idx = FindNameIndices(inputNames, requestedNames);
    matrixOut = double(matrixIn(idx, idx));
end


function T = MatrixToNamedTable(matrix, names)
    variableNames = matlab.lang.makeValidName(cellstr(names));
    T = array2table(matrix, 'VariableNames', variableNames, ...
        'RowNames', cellstr(names));
end


function T = BuildCohortTable(patientResults)
    n = numel(patientResults);
    patientID = strings(n, 1);
    nUsableParts = nan(n, 1);
    nConfiguredCorrMetabolites = nan(n, 1);
    nConfiguredCorrPairs = nan(n, 1);
    basisVsDeGraaf_matrixCorrelation = nan(n, 1);
    basisVsDeGraaf_RMSE = nan(n, 1);
    basisVsDeGraaf_MAD = nan(n, 1);
    basisVsDeGraaf_signAgreement = nan(n, 1);
    wishartDerivedRVsDeGraaf_matrixCorrelation = nan(n, 1);
    wishartDerivedRVsDeGraaf_RMSE = nan(n, 1);
    normalizedBasisVsWishart_matrixCorrelation = nan(n, 1);
    normalizedBasisVsWishart_RMSE = nan(n, 1);
    normalizedBasisVsWishart_diagonalVarianceCorrelation = nan(n, 1);
    medianConditionNumber = nan(n, 1);
    maximumConditionNumber = nan(n, 1);

    for k = 1:n
        p = patientResults(k);
        corr = p.comparisons.basisVsDeGraaf.configured;
        current = p.comparisons.currentWishartRVsDeGraaf.configured;
        covariance = p.comparisons.basisVsWishart.configured;

        patientID(k) = p.patientID;
        nUsableParts(k) = p.nUsableParts;
        nConfiguredCorrMetabolites(k) = corr.nMetabolites;
        nConfiguredCorrPairs(k) = corr.nPairs;
        basisVsDeGraaf_matrixCorrelation(k) = corr.matrixCorrelation;
        basisVsDeGraaf_RMSE(k) = corr.RMSE;
        basisVsDeGraaf_MAD(k) = corr.MAD;
        basisVsDeGraaf_signAgreement(k) = corr.signAgreement;
        wishartDerivedRVsDeGraaf_matrixCorrelation(k) = current.matrixCorrelation;
        wishartDerivedRVsDeGraaf_RMSE(k) = current.RMSE;
        normalizedBasisVsWishart_matrixCorrelation(k) = covariance.matrixCorrelation;
        normalizedBasisVsWishart_RMSE(k) = covariance.RMSE;
        normalizedBasisVsWishart_diagonalVarianceCorrelation(k) = ...
            covariance.diagonalVarianceCorrelation;
        medianConditionNumber(k) = p.conditionNumberSummary.median;
        maximumConditionNumber(k) = p.conditionNumberSummary.maximum;
    end

    T = table(patientID, nUsableParts, nConfiguredCorrMetabolites, ...
        nConfiguredCorrPairs, basisVsDeGraaf_matrixCorrelation, ...
        basisVsDeGraaf_RMSE, basisVsDeGraaf_MAD, ...
        basisVsDeGraaf_signAgreement, ...
        wishartDerivedRVsDeGraaf_matrixCorrelation, ...
        wishartDerivedRVsDeGraaf_RMSE, ...
        normalizedBasisVsWishart_matrixCorrelation, ...
        normalizedBasisVsWishart_RMSE, ...
        normalizedBasisVsWishart_diagonalVarianceCorrelation, ...
        medianConditionNumber, maximumConditionNumber);
end


function T = BuildCohortSummaryTable(cohortTable)
    excluded = ["patientID", "nUsableParts", "nConfiguredCorrMetabolites", ...
        "nConfiguredCorrPairs"];
    metricNames = setdiff(string(cohortTable.Properties.VariableNames), ...
        excluded, 'stable');
    nMetrics = numel(metricNames);
    medianValue = nan(nMetrics, 1);
    minimum = nan(nMetrics, 1);
    maximum = nan(nMetrics, 1);
    Q1 = nan(nMetrics, 1);
    Q3 = nan(nMetrics, 1);
    nPatients = zeros(nMetrics, 1);

    for k = 1:nMetrics
        x = cohortTable.(metricNames(k));
        x = x(isfinite(x));
        s = SummarizeVector(x);
        medianValue(k) = s.median;
        minimum(k) = s.minimum;
        maximum(k) = s.maximum;
        Q1(k) = s.Q1;
        Q3(k) = s.Q3;
        nPatients(k) = numel(x);
    end

    T = table(metricNames(:), nPatients, medianValue, minimum, maximum, Q1, Q3, ...
        'VariableNames', {'metric', 'nPatients', 'median', 'minimum', ...
        'maximum', 'Q1', 'Q3'});
end


function PrintDetailedPatientReport(p)
    fprintf('\nDetailed basis diagnostic: %s\n', p.patientID);
    fprintf('Parts usable: %d/%d\n', p.nUsableParts, p.nParts);
    fprintf('Spectral points: %s\n', ...
        mat2str(unique([p.partDiagnostics.spectralPointCount])));
    fprintf('Component count per part: min %d, median %.1f, max %d\n', ...
        min(p.hSummaryTable.componentCount), ...
        median(p.hSummaryTable.componentCount), ...
        max(p.hSummaryTable.componentCount));
    fprintf('Condition number: median %.6g, max %.6g\n', ...
        p.conditionNumberSummary.median, p.conditionNumberSummary.maximum);
    fprintf('All H positive definite: %d\n', all(p.hSummaryTable.positiveDefinite));
    fprintf('Configured correlation-valid names (%d): %s\n', ...
        numel(p.validCorrelationNames), strjoin(p.validCorrelationNames, ', '));
    fprintf('Configured covariance-valid names (%d): %s\n', ...
        numel(p.validCovarianceNames), strjoin(p.validCovarianceNames, ', '));

    fprintf('\nBasis vs R_DeGraaf, configured panel\n');
    disp(ComparisonSummaryTable(p.comparisons.basisVsDeGraaf.configured))
    disp(p.comparisons.basisVsDeGraaf.configured.topDiscrepancies)

    fprintf('Basis vs R_DeGraaf, final 11 panel\n');
    disp(ComparisonSummaryTable(p.comparisons.basisVsDeGraaf.final11))
    disp(p.comparisons.basisVsDeGraaf.final11.topDiscrepancies)

    fprintf('Normalized basis covariance vs current Wishart covariance\n');
    disp(CovarianceSummaryTable(p.comparisons.basisVsWishart.configured))

    fprintf('Current Wishart-derived R vs R_DeGraaf\n');
    disp(ComparisonSummaryTable( ...
        p.comparisons.currentWishartRVsDeGraaf.configured))
end


function T = ComparisonSummaryTable(c)
    T = table(c.nMetabolites, c.nPairs, c.matrixCorrelation, c.RMSE, ...
        c.MAD, c.medianAbsoluteDifference, c.meanSignedDifference, ...
        c.maximumAbsoluteDifference, c.signAgreement, ...
        'VariableNames', {'nMetabolites', 'nPairs', 'matrixCorrelation', ...
        'RMSE', 'MAD', 'medianAbsoluteDifference', 'meanSignedDifference', ...
        'maximumAbsoluteDifference', 'signAgreement'});
end


function T = CovarianceSummaryTable(c)
    T = table(c.nMetabolites, c.nEntries, c.matrixCorrelation, c.RMSE, ...
        c.MAD, c.diagonalVarianceCorrelation, c.bestFitScale, ...
        c.bestScaledRMSE, c.bestScaledMAD, ...
        'VariableNames', {'nMetabolites', 'nEntries', 'matrixCorrelation', ...
        'normalizedRMSE', 'normalizedMAD', 'diagonalVarianceCorrelation', ...
        'bestFitScale', 'bestScaledRMSE', 'bestScaledMAD'});
end


function [files, metadata] = MakeGroupDiagnosticFigures(group, opts)
    names = group.metaboliteNames;
    REmpirical = group.REmpirical;
    RDeGraaf = group.RDeGraaf;
    RBasis = group.RBasis;

    empiricalMinusDeGraaf = REmpirical - RDeGraaf;
    empiricalMinusBasis = REmpirical - RBasis;
    basisMinusDeGraaf = RBasis - RDeGraaf;
    allDifferences = [empiricalMinusDeGraaf(:); ...
        empiricalMinusBasis(:); basisMinusDeGraaf(:)];
    allDifferences = allDifferences(isfinite(allDifferences));
    if isempty(allDifferences)
        maxAbsDifference = 1;
    else
        maxAbsDifference = max(abs(allDifferences));
        if maxAbsDifference == 0
            maxAbsDifference = 1;
        end
    end
    differenceLimits = [-maxAbsDifference, maxAbsDifference];

    groupDir = fullfile(opts.outputDir, "Group");
    if opts.saveFigures && ~isfolder(groupDir)
        mkdir(groupDir);
    end

    files = strings(9, 1);
    descriptions = [ ...
        "Group empirical correlation matrix"; ...
        "Group LCModel/DeGraaf correlation matrix"; ...
        "Group basis-derived correlation matrix"; ...
        "Group empirical R minus group DeGraaf R"; ...
        "Group empirical R minus group basis R"; ...
        "Group basis R minus group DeGraaf R"; ...
        "Group DeGraaf R versus group empirical R"; ...
        "Group basis R versus group empirical R"; ...
        "Group DeGraaf R versus group basis R"];

    f = MatrixFigure(REmpirical, names, ...
        "Group empirical correlation", [-1 1], opts);
    colormap(f, gray(256));
    files(1) = SaveDiagnosticFigure(f, groupDir, ...
        "01_group_empirical_correlation.png", opts);

    f = MatrixFigure(RDeGraaf, names, ...
        "Group LCModel/DeGraaf correlation", [-1 1], opts);
    colormap(f, gray(256));
    files(2) = SaveDiagnosticFigure(f, groupDir, ...
        "02_group_DeGraaf_correlation.png", opts);

    f = MatrixFigure(RBasis, names, ...
        "Group basis-derived correlation", [-1 1], opts);
    colormap(f, gray(256));
    files(3) = SaveDiagnosticFigure(f, groupDir, ...
        "03_group_basis_correlation.png", opts);

    differenceMap = BlueWhiteRedMap(256);
    f = MatrixFigure(empiricalMinusDeGraaf, names, ...
        "Group empirical R - group DeGraaf R", differenceLimits, opts);
    colormap(f, differenceMap);
    files(4) = SaveDiagnosticFigure(f, groupDir, ...
        "04_group_empirical_minus_DeGraaf.png", opts);

    f = MatrixFigure(empiricalMinusBasis, names, ...
        "Group empirical R - group basis R", differenceLimits, opts);
    colormap(f, differenceMap);
    files(5) = SaveDiagnosticFigure(f, groupDir, ...
        "05_group_empirical_minus_basis.png", opts);

    f = MatrixFigure(basisMinusDeGraaf, names, ...
        "Group basis R - group DeGraaf R", differenceLimits, opts);
    colormap(f, differenceMap);
    files(6) = SaveDiagnosticFigure(f, groupDir, ...
        "06_group_basis_minus_DeGraaf.png", opts);

    mask = group.commonPairMask;
    empiricalValues = REmpirical(mask);
    deGraafValues = RDeGraaf(mask);
    basisValues = RBasis(mask);

    c = group.comparisons.empiricalVsDeGraaf;
    f = ScatterFigure(deGraafValues, empiricalValues, ...
        "Group DeGraaf R", "Group empirical R", ...
        "Group empirical vs DeGraaf", opts, true);
    AnnotateScatter(f, c.matrixCorrelation, c.RMSE);
    files(7) = SaveDiagnosticFigure(f, groupDir, ...
        "07_group_DeGraaf_vs_empirical_scatter.png", opts);

    c = group.comparisons.empiricalVsBasis;
    f = ScatterFigure(basisValues, empiricalValues, ...
        "Group basis R", "Group empirical R", ...
        "Group empirical vs basis", opts, true);
    AnnotateScatter(f, c.matrixCorrelation, c.RMSE);
    files(8) = SaveDiagnosticFigure(f, groupDir, ...
        "08_group_basis_vs_empirical_scatter.png", opts);

    c = group.comparisons.basisVsDeGraaf;
    f = ScatterFigure(deGraafValues, basisValues, ...
        "Group DeGraaf R", "Group basis R", ...
        "Group basis vs DeGraaf", opts, true);
    AnnotateScatter(f, c.matrixCorrelation, c.RMSE);
    files(9) = SaveDiagnosticFigure(f, groupDir, ...
        "09_group_DeGraaf_vs_basis_scatter.png", opts);

    metadata = table((1:9).', descriptions, files, ...
        'VariableNames', {'figureNumber', 'description', 'file'});
    metadata.metabolitePanel = repmat(strjoin(names, ", "), 9, 1);
    metadata.differenceColorLimit = repmat(maxAbsDifference, 9, 1);
end


function AnnotateScatter(fig, matrixCorrelation, rmse)
    ax = findobj(fig, 'Type', 'axes');
    if isempty(ax)
        return;
    end
    text(ax(1), 0.04, 0.96, ...
        sprintf('Pearson r = %.3f\nRMSE = %.3f', matrixCorrelation, rmse), ...
        'Units', 'normalized', 'VerticalAlignment', 'top', ...
        'BackgroundColor', 'w', 'EdgeColor', [0.5 0.5 0.5], ...
        'Color', 'k', 'FontSize', 11);
end


function map = BlueWhiteRedMap(nColors)
    positions = [0; 0.5; 1];
    anchors = [0.12 0.32 0.80; 1 1 1; 0.80 0.12 0.12];
    map = interp1(positions, anchors, linspace(0, 1, nColors), 'linear');
end


function PrintGroupComparisonReport(group)
    mask = triu(true(group.nMetabolites), 1);
    empiricalCounts = group.patientCountMatrices.empirical(mask);
    deGraafCounts = group.patientCountMatrices.deGraaf(mask);
    basisCounts = group.patientCountMatrices.basis(mask);

    fprintf('\nGroup correlation matrix diagnostic\n');
    fprintf('Patients available: %d\n', group.totalPatientsAvailable);
    fprintf('Panel (%d): %s\n', group.nMetabolites, ...
        strjoin(group.metaboliteNames, ', '));
    fprintf('Common finite off-diagonal pairs: %d\n', ...
        group.nCommonOffDiagonalPairs);
    fprintf('Patient counts per empirical element: %d to %d\n', ...
        min(empiricalCounts), max(empiricalCounts));
    fprintf('Patient counts per DeGraaf element: %d to %d\n', ...
        min(deGraafCounts), max(deGraafCounts));
    fprintf('Patient counts per basis element: %d to %d\n', ...
        min(basisCounts), max(basisCounts));

    labels = ["Empirical vs DeGraaf"; "Empirical vs basis"; ...
        "Basis vs DeGraaf"];
    values = {group.comparisons.empiricalVsDeGraaf; ...
        group.comparisons.empiricalVsBasis; ...
        group.comparisons.basisVsDeGraaf};
    matrixCorrelation = cellfun(@(x)x.matrixCorrelation, values);
    RMSE = cellfun(@(x)x.RMSE, values);
    MAD = cellfun(@(x)x.MAD, values);
    meanSignedDifference = cellfun(@(x)x.meanSignedDifference, values);
    T = table(labels, matrixCorrelation, RMSE, MAD, meanSignedDifference, ...
        'VariableNames', {'comparison', 'matrixCorrelation', 'RMSE', ...
        'MAD', 'meanSignedDifference'});
    disp(T)
end


function [files, metadata] = MakeDiagnosticFigures(p, opts)
    names = p.validCorrelationNames;
    idx = FindNameIndices(p.reportedNames, names);
    RBasis = p.RBasisReported(idx, idx);
    RDeGraaf = p.RDeGraaf(idx, idx);
    RWishart = p.RFromCurrentWishart(idx, idx);

    covarianceNames = p.validCovarianceNames;
    covarianceComparison = p.comparisons.basisVsWishart.configured;
    CBasisNorm = covarianceComparison.normalizedBasis;
    CWishartNorm = covarianceComparison.normalizedWishart;

    patientDir = fullfile(opts.outputDir, p.patientID);
    if opts.saveFigures && ~isfolder(patientDir)
        mkdir(patientDir);
    end

    files = strings(6, 1);
    descriptions = [ ...
        "Heatmap R_basis"; ...
        "Heatmap R_DeGraaf"; ...
        "Scatter R_DeGraaf vs R_basis"; ...
        "Heatmap correlation(C_WishartCurrent)"; ...
        "Scatter R_DeGraaf vs correlation(C_WishartCurrent)"; ...
        "Scatter normalized C_WishartCurrent vs C_basis"];

    f1 = MatrixFigure(RBasis, names, "R_basis: " + p.patientID, [-1 1], opts);
    files(1) = SaveDiagnosticFigure(f1, patientDir, "01_R_basis_heatmap.png", opts);

    f2 = MatrixFigure(RDeGraaf, names, "R_DeGraaf: " + p.patientID, [-1 1], opts);
    files(2) = SaveDiagnosticFigure(f2, patientDir, "02_R_DeGraaf_heatmap.png", opts);

    [x, y] = UpperPairVectors(RDeGraaf, RBasis, false);
    f3 = ScatterFigure(x, y, "R_DeGraaf", "R_basis", ...
        "Basis vs pairwise correlation: " + p.patientID, opts, true);
    files(3) = SaveDiagnosticFigure(f3, patientDir, "03_basis_vs_DeGraaf_scatter.png", opts);

    f4 = MatrixFigure(RWishart, names, ...
        "correlation(C_WishartCurrent): " + p.patientID, [-1 1], opts);
    files(4) = SaveDiagnosticFigure(f4, patientDir, ...
        "04_R_from_current_Wishart_heatmap.png", opts);

    [x, y] = UpperPairVectors(RDeGraaf, RWishart, false);
    f5 = ScatterFigure(x, y, "R_DeGraaf", "correlation(C_WishartCurrent)", ...
        "Two current model summaries: " + p.patientID, opts, true);
    files(5) = SaveDiagnosticFigure(f5, patientDir, ...
        "05_current_Wishart_R_vs_DeGraaf_scatter.png", opts);

    [x, y] = UpperPairVectors(CWishartNorm, CBasisNorm, true);
    f6 = ScatterFigure(x, y, "C_WishartCurrent / trace", ...
        "C_basis / trace", "Normalized covariance: " + p.patientID, opts, false);
    files(6) = SaveDiagnosticFigure(f6, patientDir, ...
        "06_normalized_covariance_scatter.png", opts);

    metadata = table((1:6).', descriptions, files, ...
        'VariableNames', {'figureNumber', 'description', 'file'});
    metadata.correlationMetabolites = repmat(strjoin(names, ", "), 6, 1);
    metadata.covarianceMetabolites = repmat(strjoin(covarianceNames, ", "), 6, 1);
end


function fig = MatrixFigure(A, names, titleText, limits, opts)
    fig = figure('Visible', char(opts.figureVisible), 'Color', 'w', ...
        'Position', [100 100 900 780]);
    imagesc(A);
    axis image;
    cb = colorbar;
    clim(limits);
    colormap(parula(256));
    title(titleText, 'Interpreter', 'none', 'Color', 'k');
    xticks(1:numel(names));
    yticks(1:numel(names));
    xticklabels(names);
    yticklabels(names);
    xtickangle(45);
    set(gca, 'TickLabelInterpreter', 'none', 'FontSize', 8, ...
        'Color', 'w', 'XColor', 'k', 'YColor', 'k');
    cb.Color = 'k';
end


function fig = ScatterFigure(x, y, xLabelText, yLabelText, titleText, opts, fixedCorrelationLimits)
    finite = isfinite(x) & isfinite(y);
    x = x(finite);
    y = y(finite);
    fig = figure('Visible', char(opts.figureVisible), 'Color', 'w', ...
        'Position', [100 100 760 680]);
    scatter(x, y, 28, 'filled', 'MarkerFaceAlpha', 0.65);
    hold on;
    if fixedCorrelationLimits
        limits = [-1 1];
    else
        allValues = [x(:); y(:)];
        if isempty(allValues)
            limits = [-1 1];
        else
            limits = [min(allValues), max(allValues)];
            if limits(1) == limits(2)
                limits = limits + [-1 1] * max(eps, abs(limits(1)) * 0.05);
            end
        end
    end
    plot(limits, limits, 'k--', 'LineWidth', 1.2);
    xlim(limits);
    ylim(limits);
    axis square;
    grid on;
    ax = gca;
    ax.Color = 'w';
    ax.XColor = 'k';
    ax.YColor = 'k';
    ax.GridColor = [0.75 0.75 0.75];
    xlabel(xLabelText, 'Interpreter', 'none', 'Color', 'k');
    ylabel(yLabelText, 'Interpreter', 'none', 'Color', 'k');
    title(titleText, 'Interpreter', 'none', 'Color', 'k');
end


function [x, y] = UpperPairVectors(A, B, includeDiagonal)
    if includeDiagonal
        mask = triu(true(size(A)));
    else
        mask = triu(true(size(A)), 1);
    end
    x = A(mask);
    y = B(mask);
end


function file = SaveDiagnosticFigure(fig, directory, filename, opts)
    if opts.saveFigures
        file = string(fullfile(directory, filename));
        exportgraphics(fig, file, 'Resolution', 200);
    else
        file = "";
    end
    if opts.closeFigures
        close(fig);
    end
end


function names = DefaultReportedNames()
    names = [ ...
        "NAA"; "NAAG"; "Cr"; "PCr"; "GPC"; "PCh"; "Glu"; "Gln"; ...
        "GABA"; "GSH"; "Tau"; "Asc"; "Glc"; "Ace"; "mI"; "sI"; ...
        "Asp"; "Lac"; "GPC+PCh"; "NAA+NAAG"; "Cr+PCr"; "Glu+Gln"];
end


function names = DefaultFinalPanel()
    names = ["GABA"; "GSH"; "Tau"; "mI"; "sI"; "Asp"; "Lac"; ...
        "GPC+PCh"; "NAA+NAAG"; "Cr+PCr"; "Glu+Gln"];
end


function names = DefaultIndependentCandidateNames()
    names = [ ...
        "NAA"; "NAAG"; "Cr"; "PCr"; "GPC"; "PCh"; "Glu"; "Gln"; ...
        "GABA"; "GSH"; "Tau"; "Asc"; "Glc"; "Ace"; "mI"; "sI"; ...
        "Asp"; "Lac"; "Lip13a"; "Lip13b"; "Lip09"; "MM09"; ...
        "Lip20"; "MM20"; "MM12"; "MM14"; "MM17"; "-CrCH2"];
end


function tf = HasFieldOrProperty(value, name)
    name = char(name);
    tf = (isstruct(value) && isfield(value, name)) || ...
        (isobject(value) && isprop(value, name));
end


function value = GetFieldOrProperty(container, name)
    value = container.(char(name));
end
