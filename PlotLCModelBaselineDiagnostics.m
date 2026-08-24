function diagnostics = PlotLCModelBaselineDiagnostics( ...
    covOutputs, baselineDiagCfg, pairwiseView)
%PlotLCModelBaselineDiagnostics Inspect Division-1 fitted LCModel baselines.
%
% This diagnostic reuses the fitData and matching .coord filenames already
% retained in covOutputs. It validates every loaded patient, keeps detailed
% figures patient-specific, and optionally summarizes patient-wise PCA.

    arguments
        covOutputs (1, 1) struct
        baselineDiagCfg (1, 1) struct
        pairwiseView (1, 1) struct = struct()
    end

    if ~isfield(baselineDiagCfg, 'patientID')
        error('baselineDiagCfg.patientID is required.');
    end
    if ~isfield(covOutputs, 'patientResults') || isempty(covOutputs.patientResults)
        error('covOutputs.patientResults is empty or unavailable.');
    end

    selectedPatientID = string(baselineDiagCfg.patientID);
    if ~isscalar(selectedPatientID) || ismissing(selectedPatientID) || ...
            strlength(selectedPatientID) == 0
        error('baselineDiagCfg.patientID must be one nonempty patient ID.');
    end

    doPCA = ReadLogicalOption(baselineDiagCfg, 'doPCA', false);
    doCohortPCA = ReadLogicalOption( ...
        baselineDiagCfg, 'doCohortPCA', false);
    nPCsToPlot = ReadPositiveIntegerOption( ...
        baselineDiagCfg, 'nPCsToPlot', 3);
    doMetaboliteRegression = ReadLogicalOption( ...
        baselineDiagCfg, 'doMetaboliteRegression', false);
    regressionMinObservations = ReadPositiveIntegerOption( ...
        baselineDiagCfg, 'regressionMinObservations', 8);
    if regressionMinObservations <= 4
        error(['baselineDiagCfg.regressionMinObservations must exceed four ', ...
            'for three predictors plus an intercept.']);
    end
    regressionIgnoreZeros = ReadLogicalOption( ...
        baselineDiagCfg, 'regressionIgnoreZeros', true);
    regressionExamplePatientID = ReadStringOption( ...
        baselineDiagCfg, 'regressionExamplePatientID', selectedPatientID);
    regressionExampleMetabolite = ReadStringOption( ...
        baselineDiagCfg, 'regressionExampleMetabolite', "NAA");
    doPairResidualization = ReadLogicalOption( ...
        baselineDiagCfg, 'doPairResidualization', false);
    pairExamplePatientID = ReadStringOption( ...
        baselineDiagCfg, 'pairExamplePatientID', selectedPatientID);
    pairExampleMetabolites = ReadStringPairOption( ...
        baselineDiagCfg, 'pairExampleMetabolites', ...
        ["NAA+NAAG", "Cr+PCr"]);
    if doPairResidualization && ~doMetaboliteRegression
        error(['Pair residualization requires ', ...
            'baselineDiagCfg.doMetaboliteRegression = true.']);
    end
    if doPairResidualization
        ValidatePairwiseView(pairwiseView);
    end

    selectedBaselineMatrix = [];
    selectedPPM = [];
    selectedPCA = struct();

    nPatients = numel(covOutputs.patientResults);
    patientIDs = strings(nPatients, 1);
    PC1ExplainedPct = nan(nPatients, 1);
    PC2ExplainedPct = nan(nPatients, 1);
    PC3ExplainedPct = nan(nPatients, 1);
    PC1CumulativePct = nan(nPatients, 1);
    PC1To2CumulativePct = nan(nPatients, 1);
    PC1To3CumulativePct = nan(nPatients, 1);
    maxPC1EnergyFraction = nan(nPatients, 1);
    partWithLargestAbsPC1Score = nan(nPatients, 1);
    meanBaselineDeviationRMS = nan(nPatients, 1);
    maxBaselineDeviationRMS = nan(nPatients, 1);
    partWithMaxBaselineDeviationRMS = nan(nPatients, 1);

    metaboliteList = strings(0, 1);
    if doMetaboliteRegression
        if ~isfield(covOutputs, 'metabList') || isempty(covOutputs.metabList)
            error('covOutputs.metabList is required for metabolite regression.');
        end
        metaboliteList = string(covOutputs.metabList(:));
    end
    nMetabolites = numel(metaboliteList);
    nRegressionRows = nPatients * nMetabolites;
    regressionPatientID = strings(nRegressionRows, 1);
    regressionMetabolite = strings(nRegressionRows, 1);
    regressionNObservations = zeros(nRegressionRows, 1);
    regressionBeta0 = nan(nRegressionRows, 1);
    regressionBetaPC1 = nan(nRegressionRows, 1);
    regressionBetaPC2 = nan(nRegressionRows, 1);
    regressionBetaPC3 = nan(nRegressionRows, 1);
    regressionR2 = nan(nRegressionRows, 1);
    regressionAdjustedR2 = nan(nRegressionRows, 1);
    regressionModelPValue = nan(nRegressionRows, 1);
    regressionStatus = strings(nRegressionRows, 1);
    regressionExample = struct('wasFound', false);
    emptySignals = struct( ...
        'patientID', "", ...
        'metaboliteList', strings(0, 1), ...
        'partNumbers', zeros(0, 1), ...
        'concentration', zeros(0, 0), ...
        'predicted', zeros(0, 0), ...
        'baselineAssociatedVariation', zeros(0, 0), ...
        'residual', zeros(0, 0), ...
        'regressionStatus', strings(0, 1));
    regressionSignals = repmat(emptySignals, nPatients, 1);

    for patientIndex = 1:numel(covOutputs.patientResults)
        patientResult = covOutputs.patientResults(patientIndex);

        if doMetaboliteRegression
            requiredFields = ["patientID", "coordFiles", "fitData", "partTable"];
        else
            requiredFields = ["patientID", "coordFiles", "fitData"];
        end
        if ~all(isfield(patientResult, cellstr(requiredFields)))
            error('Patient result %d lacks patientID, coordFiles, or fitData.', ...
                patientIndex);
        end

        patientID = string(patientResult.patientID);
        patientIDs(patientIndex) = patientID;
        [baselineMatrix, ppmAxis] = ExtractDivision1Baselines( ...
            patientID, patientResult.coordFiles, patientResult.fitData);

        isSelectedPatient = patientID == selectedPatientID;
        if doCohortPCA || doMetaboliteRegression || ...
                (doPCA && isSelectedPatient)
            patientPCA = ComputePatientBaselinePCA( ...
                baselineMatrix, ppmAxis, patientID);
        end

        if doCohortPCA
            PC1ExplainedPct(patientIndex) = patientPCA.pcaExplained(1);
            PC2ExplainedPct(patientIndex) = patientPCA.pcaExplained(2);
            PC3ExplainedPct(patientIndex) = patientPCA.pcaExplained(3);
            PC1CumulativePct(patientIndex) = patientPCA.pcaExplained(1);
            PC1To2CumulativePct(patientIndex) = ...
                sum(patientPCA.pcaExplained(1:2));
            PC1To3CumulativePct(patientIndex) = ...
                sum(patientPCA.pcaExplained(1:3));
            maxPC1EnergyFraction(patientIndex) = ...
                patientPCA.maxPC1EnergyFraction;
            partWithLargestAbsPC1Score(patientIndex) = ...
                patientPCA.partWithLargestAbsPC1Score;
            meanBaselineDeviationRMS(patientIndex) = ...
                mean(patientPCA.baselineDeviationRMS);
            maxBaselineDeviationRMS(patientIndex) = ...
                max(patientPCA.baselineDeviationRMS);
            partWithMaxBaselineDeviationRMS(patientIndex) = ...
                patientPCA.partWithMaxBaselineDeviationRMS;
        end

        if doMetaboliteRegression
            concentrationMatrix = AlignConcentrationsToParts( ...
                patientResult.partTable, metaboliteList, ...
                patientPCA.partNumbers, patientID);
            predictedMatrix = nan(36, nMetabolites);
            baselineAssociatedMatrix = nan(36, nMetabolites);
            residualMatrix = nan(36, nMetabolites);
            patientRegressionStatus = strings(nMetabolites, 1);

            for metaboliteIndex = 1:nMetabolites
                regressionIndex = ...
                    (patientIndex - 1) * nMetabolites + metaboliteIndex;
                metabolite = metaboliteList(metaboliteIndex);
                regression = FitBaselinePCMetaboliteRegression( ...
                    concentrationMatrix(:, metaboliteIndex), ...
                    patientPCA.pcaScore(:, 1:3), ...
                    patientPCA.partNumbers, ...
                    regressionMinObservations, ...
                    regressionIgnoreZeros);

                regressionPatientID(regressionIndex) = patientID;
                regressionMetabolite(regressionIndex) = metabolite;
                regressionNObservations(regressionIndex) = ...
                    regression.nObservations;
                regressionBeta0(regressionIndex) = regression.beta(1);
                regressionBetaPC1(regressionIndex) = regression.beta(2);
                regressionBetaPC2(regressionIndex) = regression.beta(3);
                regressionBetaPC3(regressionIndex) = regression.beta(4);
                regressionR2(regressionIndex) = regression.R2;
                regressionAdjustedR2(regressionIndex) = regression.adjustedR2;
                regressionModelPValue(regressionIndex) = regression.modelPValue;
                regressionStatus(regressionIndex) = regression.status;
                predictedMatrix(:, metaboliteIndex) = regression.predicted;
                baselineAssociatedMatrix(:, metaboliteIndex) = ...
                    regression.baselineAssociatedVariation;
                residualMatrix(:, metaboliteIndex) = regression.residual;
                patientRegressionStatus(metaboliteIndex) = regression.status;

                if patientID == regressionExamplePatientID && ...
                        metabolite == regressionExampleMetabolite
                    regressionExample = regression;
                    regressionExample.wasFound = true;
                    regressionExample.patientID = patientID;
                    regressionExample.metabolite = metabolite;
                end
            end

            regressionSignals(patientIndex).patientID = patientID;
            regressionSignals(patientIndex).metaboliteList = metaboliteList;
            regressionSignals(patientIndex).partNumbers = ...
                patientPCA.partNumbers;
            regressionSignals(patientIndex).concentration = ...
                concentrationMatrix;
            regressionSignals(patientIndex).predicted = predictedMatrix;
            regressionSignals(patientIndex).baselineAssociatedVariation = ...
                baselineAssociatedMatrix;
            regressionSignals(patientIndex).residual = residualMatrix;
            regressionSignals(patientIndex).regressionStatus = ...
                patientRegressionStatus;
        end

        if isSelectedPatient
            selectedBaselineMatrix = baselineMatrix;
            selectedPPM = ppmAxis;
            if doPCA || doCohortPCA || doMetaboliteRegression
                selectedPCA = patientPCA;
            end
        end
    end

    if isempty(selectedBaselineMatrix)
        error('Diagnostic patient %s was not found in covOutputs.', selectedPatientID);
    end

    meanBaseline = mean(selectedBaselineMatrix, 1);
    baselineDeviation = selectedBaselineMatrix - meanBaseline;
    pointwiseSD = std(selectedBaselineMatrix, 0, 1);

    fprintf('\nLCModel fitted-baseline diagnostic\n');
    fprintf('Patient: %s\n', selectedPatientID);
    fprintf('Number of Division-1 baselines: %d\n', size(selectedBaselineMatrix, 1));
    fprintf('Spectral points per baseline: %d\n', size(selectedBaselineMatrix, 2));
    fprintf('PPM range: %.1f to %.1f\n', selectedPPM(1), selectedPPM(end));
    fprintf('Mean pointwise SD across ppm: %.6g\n', mean(pointwiseSD));
    fprintf('Maximum pointwise SD across ppm: %.6g\n', max(pointwiseSD));

    lineColors = lines(size(selectedBaselineMatrix, 1));

    figure('Name', sprintf('%s Division-1 fitted baselines', selectedPatientID), ...
        'NumberTitle', 'off');
    colororder(lineColors);
    plot(selectedPPM, selectedBaselineMatrix.', 'LineWidth', 0.8);
    set(gca, 'XDir', 'reverse');
    xlabel('ppm');
    ylabel('Fitted baseline amplitude');
    title(sprintf('%s: all 36 Division-1 fitted baselines', selectedPatientID), ...
        'Interpreter', 'none');
    grid on;
    box on;

    figure('Name', sprintf('%s Division-1 baseline deviations', selectedPatientID), ...
        'NumberTitle', 'off');
    colororder(lineColors);
    plot(selectedPPM, baselineDeviation.', 'LineWidth', 0.8);
    set(gca, 'XDir', 'reverse');
    xlabel('ppm');
    ylabel('Baseline amplitude minus patient mean');
    title(sprintf('%s: fitted-baseline deviations from patient mean', ...
        selectedPatientID), 'Interpreter', 'none');
    grid on;
    box on;

    if doPCA
    pcaCoeff = selectedPCA.pcaCoeff;
    pcaScore = selectedPCA.pcaScore;
    pcaExplained = selectedPCA.pcaExplained;
    baselineDeviationRMS = selectedPCA.baselineDeviationRMS;
    partNumbers = selectedPCA.partNumbers;
    nAvailablePCs = numel(pcaExplained);
    if nPCsToPlot > nAvailablePCs
        error('Requested %d PCs, but only %d are available.', ...
            nPCsToPlot, nAvailablePCs);
    end

    fprintf('\nPCA explained variance:\n');
    fprintf('  PC1: %.2f %%\n', pcaExplained(1));
    fprintf('  PC2: %.2f %%\n', pcaExplained(2));
    fprintf('  PC3: %.2f %%\n', pcaExplained(3));
    fprintf('  cumulative PC1: %.2f %%\n', pcaExplained(1));
    fprintf('  cumulative PC1-2: %.2f %%\n', sum(pcaExplained(1:2)));
    fprintf('  cumulative PC1-3: %.2f %%\n', sum(pcaExplained(1:3)));
    fprintf('Maximum baseline-deviation RMS:\n');
    fprintf('  part: %d\n', selectedPCA.partWithMaxBaselineDeviationRMS);
    fprintf('  RMS: %.6g\n', max(baselineDeviationRMS));

    nExplainedToPlot = min(10, nAvailablePCs);
    figure('Name', sprintf('%s baseline PCA explained variance', ...
        selectedPatientID), 'NumberTitle', 'off');
    bar(1:nExplainedToPlot, pcaExplained(1:nExplainedToPlot));
    xlabel('Principal component');
    ylabel('Explained variance (%)');
    title(sprintf('%s: fitted-baseline PCA explained variance', ...
        selectedPatientID), 'Interpreter', 'none');
    xticks(1:nExplainedToPlot);
    grid on;
    box on;

    figure('Name', sprintf('%s baseline PCA spectral loadings', ...
        selectedPatientID), 'NumberTitle', 'off');
    loadingLayout = tiledlayout(nPCsToPlot, 1, ...
        'TileSpacing', 'compact', 'Padding', 'compact');
    for pcIndex = 1:nPCsToPlot
        nexttile;
        plot(selectedPPM, pcaCoeff(:, pcIndex), 'LineWidth', 1.1);
        set(gca, 'XDir', 'reverse');
        ylabel(sprintf('PC%d', pcIndex));
        title(sprintf('PC%d loading (%.2f%% explained)', ...
            pcIndex, pcaExplained(pcIndex)));
        grid on;
        box on;
    end
    xlabel(loadingLayout, 'ppm');
    title(loadingLayout, sprintf('%s: baseline PCA spectral loadings', ...
        selectedPatientID), 'Interpreter', 'none');

    figure('Name', sprintf('%s baseline PCA scores by part', ...
        selectedPatientID), 'NumberTitle', 'off');
    scoreLayout = tiledlayout(nPCsToPlot, 1, ...
        'TileSpacing', 'compact', 'Padding', 'compact');
    for pcIndex = 1:nPCsToPlot
        nexttile;
        plot(partNumbers, pcaScore(:, pcIndex), '-o', ...
            'LineWidth', 1.0, 'MarkerSize', 3);
        xlim([1, 36]);
        xticks(1:5:36);
        ylabel(sprintf('PC%d score', pcIndex));
        title(sprintf('PC%d scores (%.2f%% explained)', ...
            pcIndex, pcaExplained(pcIndex)));
        grid on;
        box on;
    end
    xlabel(scoreLayout, 'Division-1 part number');
    title(scoreLayout, sprintf('%s: baseline PCA scores across parts', ...
        selectedPatientID), 'Interpreter', 'none');

    figure('Name', sprintf('%s baseline-deviation RMS by part', ...
        selectedPatientID), 'NumberTitle', 'off');
    plot(partNumbers, baselineDeviationRMS, '-o', ...
        'LineWidth', 1.1, 'MarkerSize', 4);
    xlim([1, 36]);
    xticks(1:5:36);
    xlabel('Division-1 part number');
    ylabel('Baseline-deviation RMS');
    title(sprintf('%s: baseline-deviation magnitude by part', ...
        selectedPatientID), 'Interpreter', 'none');
    grid on;
    box on;
    end

    cohortTable = table();
    if doCohortPCA
        cohortTable = table( ...
            patientIDs, ...
            PC1ExplainedPct, ...
            PC2ExplainedPct, ...
            PC3ExplainedPct, ...
            PC1CumulativePct, ...
            PC1To2CumulativePct, ...
            PC1To3CumulativePct, ...
            maxPC1EnergyFraction, ...
            partWithLargestAbsPC1Score, ...
            meanBaselineDeviationRMS, ...
            maxBaselineDeviationRMS, ...
            partWithMaxBaselineDeviationRMS, ...
            'VariableNames', { ...
            'patientID', ...
            'PC1ExplainedPct', ...
            'PC2ExplainedPct', ...
            'PC3ExplainedPct', ...
            'PC1CumulativePct', ...
            'PC1To2CumulativePct', ...
            'PC1To3CumulativePct', ...
            'maxPC1EnergyFraction', ...
            'partWithLargestAbsPC1Score', ...
            'meanBaselineDeviationRMS', ...
            'maxBaselineDeviationRMS', ...
            'partWithMaxBaselineDeviationRMS'});

        PrintCohortPCASummary(cohortTable);
        PlotCohortPCAFigures(cohortTable);
    end

    metaboliteRegressionTable = table();
    metaboliteRegressionSummaryTable = table();
    regressionOverallSummary = struct();
    if doMetaboliteRegression
        metaboliteRegressionTable = table( ...
            regressionPatientID, ...
            regressionMetabolite, ...
            regressionNObservations, ...
            regressionBeta0, ...
            regressionBetaPC1, ...
            regressionBetaPC2, ...
            regressionBetaPC3, ...
            regressionR2, ...
            regressionAdjustedR2, ...
            regressionModelPValue, ...
            regressionStatus, ...
            'VariableNames', { ...
            'patientID', ...
            'metabolite', ...
            'nObservations', ...
            'beta0', ...
            'betaPC1', ...
            'betaPC2', ...
            'betaPC3', ...
            'R2', ...
            'adjustedR2', ...
            'modelPValue', ...
            'status'});

        [metaboliteRegressionSummaryTable, regressionOverallSummary] = ...
            SummarizeMetaboliteRegressions( ...
            metaboliteRegressionTable, metaboliteList);
        PrintMetaboliteRegressionSummary( ...
            metaboliteRegressionSummaryTable, regressionOverallSummary, ...
            nPatients);

        if ~regressionExample.wasFound
            error('Regression example %s / %s was not found.', ...
                regressionExamplePatientID, regressionExampleMetabolite);
        end
        PlotMetaboliteRegressionFigures( ...
            metaboliteRegressionTable, metaboliteList, regressionExample);
    end

    pairResidualizationTable = table();
    pairResidualizationSummaryTable = table();
    pairResidualizationGlobalSummary = struct();
    pairResidualizationExample = struct();
    if doPairResidualization
        [pairResidualizationTable, pairResidualizationSummaryTable, ...
            pairResidualizationGlobalSummary, pairResidualizationExample] = ...
            BuildPairResidualizationDiagnostics( ...
            pairwiseView, regressionSignals, pairExamplePatientID, ...
            pairExampleMetabolites);
        PrintPairResidualizationSummary( ...
            pairResidualizationSummaryTable, ...
            pairResidualizationGlobalSummary);
        PlotPairResidualizationFigures( ...
            pairResidualizationSummaryTable, ...
            pairResidualizationExample);
    end

    diagnostics = struct();
    diagnostics.patientTable = cohortTable;
    diagnostics.metaboliteRegressionTable = metaboliteRegressionTable;
    diagnostics.metaboliteRegressionSummaryTable = ...
        metaboliteRegressionSummaryTable;
    diagnostics.regressionOverallSummary = regressionOverallSummary;
    diagnostics.regressionExample = regressionExample;
    diagnostics.pairResidualizationTable = pairResidualizationTable;
    diagnostics.pairResidualizationSummaryTable = ...
        pairResidualizationSummaryTable;
    diagnostics.pairResidualizationGlobalSummary = ...
        pairResidualizationGlobalSummary;
    diagnostics.pairResidualizationExample = pairResidualizationExample;
    diagnostics.selectedPatient.patientID = selectedPatientID;
    diagnostics.selectedPatient.baselineMatrix = selectedBaselineMatrix;
    diagnostics.selectedPatient.meanBaseline = meanBaseline;
    diagnostics.selectedPatient.baselineDeviation = baselineDeviation;
    diagnostics.selectedPatient.ppm = selectedPPM;
    diagnostics.selectedPatient.partNumbers = (1:36).';
    if doPCA || doCohortPCA || doMetaboliteRegression
        diagnostics.selectedPatient.pcaCoeff = selectedPCA.pcaCoeff;
        diagnostics.selectedPatient.pcaScore = selectedPCA.pcaScore;
        diagnostics.selectedPatient.pcaExplained = selectedPCA.pcaExplained;
        diagnostics.selectedPatient.baselineDeviationRMS = ...
            selectedPCA.baselineDeviationRMS;
    end
end


function value = ReadLogicalOption(options, fieldName, defaultValue)
    value = defaultValue;
    if isfield(options, fieldName)
        value = options.(fieldName);
    end
    if ~islogical(value) || ~isscalar(value)
        error('baselineDiagCfg.%s must be a logical scalar.', fieldName);
    end
end


function value = ReadPositiveIntegerOption(options, fieldName, defaultValue)
    value = defaultValue;
    if isfield(options, fieldName)
        value = options.(fieldName);
    end
    if ~isnumeric(value) || ~isscalar(value) || ~isfinite(value) || ...
            value < 1 || value ~= fix(value)
        error('baselineDiagCfg.%s must be a positive integer.', fieldName);
    end
end


function value = ReadStringOption(options, fieldName, defaultValue)
    value = string(defaultValue);
    if isfield(options, fieldName)
        value = string(options.(fieldName));
    end
    if ~isscalar(value) || ismissing(value) || strlength(value) == 0
        error('baselineDiagCfg.%s must be one nonempty string.', fieldName);
    end
end


function value = ReadStringPairOption(options, fieldName, defaultValue)
    value = string(defaultValue);
    if isfield(options, fieldName)
        value = string(options.(fieldName));
    end
    value = reshape(value, 1, []);
    if numel(value) ~= 2 || any(ismissing(value)) || any(strlength(value) == 0) || ...
            value(1) == value(2)
        error('baselineDiagCfg.%s must contain two distinct metabolite names.', ...
            fieldName);
    end
end


function ValidatePairwiseView(pairwiseView)
    requiredFields = [ ...
        "patientEligibilityTable", ...
        "pairTable", ...
        "minValidParts", ...
        "minPatientsForGroupTest"];
    if ~all(isfield(pairwiseView, cellstr(requiredFields)))
        error(['Pair residualization requires the filtered analysisData.pairwise ', ...
            'view with eligibility tables and minimum-count rules.']);
    end
    requiredColumns = [ ...
        "patientID", "metaboliteA", "metaboliteB", ...
        "nValidParts", "eligible", "rEmpirical", "rModel", "reason"];
    if ~istable(pairwiseView.patientEligibilityTable) || ...
            ~all(ismember(requiredColumns, string( ...
            pairwiseView.patientEligibilityTable.Properties.VariableNames)))
        error('analysisData.pairwise.patientEligibilityTable has an invalid schema.');
    end
end


function result = ComputePatientBaselinePCA(baselineMatrix, ppmAxis, patientID)
    if size(baselineMatrix, 1) ~= 36
        error('Patient %s PCA requires exactly 36 baseline rows.', patientID);
    end

    meanBaseline = mean(baselineMatrix, 1);
    baselineDeviation = baselineMatrix - meanBaseline;

    % Observations are parts 1:36 and variables are ppm positions. X is
    % mean-centered only; no ppm-wise scaling or standardization is applied.
    [pcaCoeff, pcaScore, pcaLatent] = pca( ...
        baselineDeviation, 'Centered', false, 'Economy', true);

    totalVariance = sum(pcaLatent);
    if ~isfinite(totalVariance) || totalVariance <= 0
        error('Patient %s has no finite positive baseline variance for PCA.', ...
            patientID);
    end
    pcaExplained = 100 .* pcaLatent ./ totalVariance;

    if numel(pcaExplained) < 3
        error('Patient %s PCA returned fewer than three components.', patientID);
    end
    if size(pcaCoeff, 1) ~= numel(ppmAxis)
        error('Patient %s PCA loading length does not match its ppm grid.', ...
            patientID);
    end
    if size(pcaScore, 1) ~= 36
        error('Patient %s PCA scores do not preserve all parts 1:36.', ...
            patientID);
    end
    if any(~isfinite(pcaExplained)) || ...
            abs(sum(pcaExplained) - 100) > 1e-8
        error(['Patient %s PCA explained variance is invalid or does not ', ...
            'sum to 100%%.'], patientID);
    end

    partNumbers = (1:36).';
    pc1Score = pcaScore(:, 1);
    pc1Energy = sum(pc1Score.^2);
    if ~isfinite(pc1Energy) || pc1Energy <= 0
        error('Patient %s has invalid zero PC1 score energy.', patientID);
    end
    [~, largestPC1Index] = max(abs(pc1Score));

    baselineDeviationRMS = sqrt(mean(baselineDeviation.^2, 2));
    [~, largestRMSIndex] = max(baselineDeviationRMS);

    result = struct();
    result.meanBaseline = meanBaseline;
    result.baselineDeviation = baselineDeviation;
    result.pcaCoeff = pcaCoeff;
    result.pcaScore = pcaScore;
    result.pcaExplained = pcaExplained;
    result.partNumbers = partNumbers;
    result.maxPC1EnergyFraction = max(pc1Score.^2) ./ pc1Energy;
    result.partWithLargestAbsPC1Score = partNumbers(largestPC1Index);
    result.baselineDeviationRMS = baselineDeviationRMS;
    result.partWithMaxBaselineDeviationRMS = partNumbers(largestRMSIndex);
end


function concentrationMatrix = AlignConcentrationsToParts( ...
    partTable, metaboliteList, targetPartNumbers, patientID)

    if ~istable(partTable) || isempty(partTable) || ...
            ~ismember('part', partTable.Properties.VariableNames)
        error('Patient %s has no usable empirical partTable.', patientID);
    end

    partNumbers = double(partTable.part(:));
    if any(~isfinite(partNumbers)) || numel(unique(partNumbers)) ~= numel(partNumbers)
        error('Patient %s empirical part numbers are nonfinite or duplicated.', ...
            patientID);
    end

    [partsMatched, rowLocations] = ismember( ...
        double(targetPartNumbers(:)), partNumbers);
    if ~all(partsMatched)
        error('Patient %s empirical concentrations do not contain parts 1:36.', ...
            patientID);
    end

    variableNames = string(partTable.Properties.VariableNames);
    concentrationMatrix = nan(numel(targetPartNumbers), numel(metaboliteList));

    for metaboliteIndex = 1:numel(metaboliteList)
        columnName = string(matlab.lang.makeValidName( ...
            char(metaboliteList(metaboliteIndex))));
        columnIndex = find(variableNames == columnName, 1, 'first');
        if isempty(columnIndex)
            error('Patient %s is missing concentration column %s.', ...
                patientID, metaboliteList(metaboliteIndex));
        end
        if ~isnumeric(partTable{:, columnIndex})
            error('Patient %s concentration column %s is nonnumeric.', ...
                patientID, metaboliteList(metaboliteIndex));
        end

        values = double(partTable{:, columnIndex});
        concentrationMatrix(:, metaboliteIndex) = values(rowLocations);
    end
end


function result = FitBaselinePCMetaboliteRegression( ...
    concentrations, pcaScores, partNumbers, minObservations, ignoreZeros)

    concentrations = double(concentrations(:));
    pcaScores = double(pcaScores);
    partNumbers = double(partNumbers(:));

    if numel(concentrations) ~= 36 || ...
            ~isequal(size(pcaScores), [36, 3]) || ...
            ~isequal(partNumbers, (1:36).')
        error('Regression inputs are not aligned to Division-1 parts 1:36.');
    end

    valid = isfinite(concentrations) & all(isfinite(pcaScores), 2);
    if ignoreZeros
        valid = valid & concentrations ~= 0;
    end

    result = struct();
    result.nObservations = sum(valid);
    result.beta = nan(4, 1);
    result.R2 = NaN;
    result.adjustedR2 = NaN;
    result.modelPValue = NaN;
    result.status = "UNAVAILABLE";
    result.partNumbers = partNumbers;
    result.observed = concentrations;
    result.observed(~valid) = NaN;
    result.predicted = nan(36, 1);
    result.baselineAssociatedVariation = nan(36, 1);
    result.residual = nan(36, 1);
    result.validRows = valid;

    if result.nObservations < minObservations
        result.status = "INSUFFICIENT_OBSERVATIONS";
        return;
    end

    design = [ones(result.nObservations, 1), pcaScores(valid, :)];
    response = concentrations(valid);
    if rank(design) < 4
        result.status = "RANK_DEFICIENT";
        return;
    end

    responseSS = sum((response - mean(response)).^2);
    if ~isfinite(responseSS) || responseSS <= 0
        result.status = "CONSTANT_RESPONSE";
        return;
    end

    beta = design \ response;
    fitted = design * beta;
    residual = response - fitted;
    errorSS = sum(residual.^2);
    R2 = 1 - errorSS ./ responseSS;
    nPredictors = 3;
    errorDF = result.nObservations - nPredictors - 1;
    adjustedR2 = 1 - (1 - R2) .* ...
        (result.nObservations - 1) ./ errorDF;

    modelSS = max(responseSS - errorSS, 0);
    if errorSS <= eps(responseSS)
        modelPValue = 0;
    else
        F = (modelSS ./ nPredictors) ./ (errorSS ./ errorDF);
        modelPValue = fcdf(F, nPredictors, errorDF, 'upper');
    end

    result.beta = beta;
    result.R2 = R2;
    result.adjustedR2 = adjustedR2;
    result.modelPValue = modelPValue;
    result.status = "OK";
    result.predicted = [ones(36, 1), pcaScores] * beta;
    result.baselineAssociatedVariation = pcaScores * beta(2:4);
    result.baselineAssociatedVariation(~valid) = NaN;
    result.residual(valid) = response - fitted;
end


function [summaryTable, overallSummary] = SummarizeMetaboliteRegressions( ...
    regressionTable, metaboliteList)

    validRows = regressionTable.status == "OK" & ...
        isfinite(regressionTable.adjustedR2);
    validValues = regressionTable.adjustedR2(validRows);
    if isempty(validValues)
        error('No valid baseline-PC metabolite regressions were available.');
    end

    nMetabolites = numel(metaboliteList);
    nPatients = zeros(nMetabolites, 1);
    medianAdjustedR2 = nan(nMetabolites, 1);
    meanAdjustedR2 = nan(nMetabolites, 1);
    Q1AdjustedR2 = nan(nMetabolites, 1);
    Q3AdjustedR2 = nan(nMetabolites, 1);
    minAdjustedR2 = nan(nMetabolites, 1);
    maxAdjustedR2 = nan(nMetabolites, 1);
    nAdjustedR2AtLeast010 = zeros(nMetabolites, 1);
    nAdjustedR2AtLeast025 = zeros(nMetabolites, 1);
    nAdjustedR2AtLeast050 = zeros(nMetabolites, 1);

    for metaboliteIndex = 1:nMetabolites
        rows = validRows & ...
            regressionTable.metabolite == metaboliteList(metaboliteIndex);
        values = regressionTable.adjustedR2(rows);
        if isempty(values)
            continue;
        end
        quartiles = prctile(values, [25, 75]);
        nPatients(metaboliteIndex) = numel(values);
        medianAdjustedR2(metaboliteIndex) = median(values);
        meanAdjustedR2(metaboliteIndex) = mean(values);
        Q1AdjustedR2(metaboliteIndex) = quartiles(1);
        Q3AdjustedR2(metaboliteIndex) = quartiles(2);
        minAdjustedR2(metaboliteIndex) = min(values);
        maxAdjustedR2(metaboliteIndex) = max(values);
        nAdjustedR2AtLeast010(metaboliteIndex) = sum(values >= 0.10);
        nAdjustedR2AtLeast025(metaboliteIndex) = sum(values >= 0.25);
        nAdjustedR2AtLeast050(metaboliteIndex) = sum(values >= 0.50);
    end

    summaryTable = table( ...
        metaboliteList, ...
        nPatients, ...
        medianAdjustedR2, ...
        meanAdjustedR2, ...
        Q1AdjustedR2, ...
        Q3AdjustedR2, ...
        minAdjustedR2, ...
        maxAdjustedR2, ...
        nAdjustedR2AtLeast010, ...
        nAdjustedR2AtLeast025, ...
        nAdjustedR2AtLeast050, ...
        'VariableNames', { ...
        'metabolite', ...
        'nPatients', ...
        'medianAdjustedR2', ...
        'meanAdjustedR2', ...
        'Q1AdjustedR2', ...
        'Q3AdjustedR2', ...
        'minAdjustedR2', ...
        'maxAdjustedR2', ...
        'nAdjustedR2AtLeast010', ...
        'nAdjustedR2AtLeast025', ...
        'nAdjustedR2AtLeast050'});

    overallQuartiles = prctile(validValues, [25, 75]);
    overallSummary = struct();
    overallSummary.totalRegressions = height(regressionTable);
    overallSummary.nValidRegressions = numel(validValues);
    overallSummary.medianAdjustedR2 = median(validValues);
    overallSummary.meanAdjustedR2 = mean(validValues);
    overallSummary.Q1AdjustedR2 = overallQuartiles(1);
    overallSummary.Q3AdjustedR2 = overallQuartiles(2);
end


function PrintMetaboliteRegressionSummary( ...
    summaryTable, overallSummary, nPatients)

    fprintf('\nBaseline-PC metabolite regression\n');
    fprintf('Patients: %d\n', nPatients);
    fprintf('PCs used: 3\n');
    fprintf('Valid regressions: %d of %d\n', ...
        overallSummary.nValidRegressions, overallSummary.totalRegressions);
    fprintf('Overall adjusted R-squared:\n');
    fprintf('  median: %.4f\n', overallSummary.medianAdjustedR2);
    fprintf('  mean: %.4f\n', overallSummary.meanAdjustedR2);
    fprintf('  Q1: %.4f\n', overallSummary.Q1AdjustedR2);
    fprintf('  Q3: %.4f\n', overallSummary.Q3AdjustedR2);

    [~, order] = sort(summaryTable.medianAdjustedR2, 'descend', ...
        'MissingPlacement', 'last');
    fprintf('Metabolite summary by median adjusted R-squared:\n');
    fprintf('  metabolite  n  median  mean  Q1  Q3  >=.10  >=.25  >=.50\n');
    for rowIndex = reshape(order, 1, [])
        fprintf('  %-11s %2d  %.4f  %.4f  %.4f  %.4f  %3d  %3d  %3d\n', ...
            summaryTable.metabolite(rowIndex), ...
            summaryTable.nPatients(rowIndex), ...
            summaryTable.medianAdjustedR2(rowIndex), ...
            summaryTable.meanAdjustedR2(rowIndex), ...
            summaryTable.Q1AdjustedR2(rowIndex), ...
            summaryTable.Q3AdjustedR2(rowIndex), ...
            summaryTable.nAdjustedR2AtLeast010(rowIndex), ...
            summaryTable.nAdjustedR2AtLeast025(rowIndex), ...
            summaryTable.nAdjustedR2AtLeast050(rowIndex));
    end
end


function PlotMetaboliteRegressionFigures( ...
    regressionTable, metaboliteList, regressionExample)

    validRows = regressionTable.status == "OK" & ...
        isfinite(regressionTable.adjustedR2);
    nValid = sum(validRows);
    adjustedR2 = nan(nValid, 1);
    groupIndex = nan(nValid, 1);
    nextRow = 1;
    for metaboliteIndex = 1:numel(metaboliteList)
        rows = find(validRows & ...
            regressionTable.metabolite == metaboliteList(metaboliteIndex));
        rowRange = nextRow:(nextRow + numel(rows) - 1);
        adjustedR2(rowRange) = regressionTable.adjustedR2(rows);
        groupIndex(rowRange) = metaboliteIndex;
        nextRow = nextRow + numel(rows);
    end

    figure('Name', 'Baseline-PC adjusted R-squared by metabolite', ...
        'NumberTitle', 'off');
    metaboliteGroups = categorical( ...
        groupIndex, 1:numel(metaboliteList), cellstr(metaboliteList));
    boxchart(metaboliteGroups, adjustedR2, 'MarkerStyle', 'none');
    yline(0, '--');
    xtickangle(45);
    set(gca, 'TickLabelInterpreter', 'none');
    xlabel('Metabolite');
    ylabel('Adjusted R-squared');
    title('Association of metabolite concentrations with baseline PC1-PC3');
    grid on;
    box on;

    figure('Name', sprintf('%s %s baseline-PC regression example', ...
        regressionExample.patientID, regressionExample.metabolite), ...
        'NumberTitle', 'off');
    plot(regressionExample.partNumbers, regressionExample.observed, '-o', ...
        'LineWidth', 1.1, 'MarkerSize', 4);
    hold on;
    plot(regressionExample.partNumbers, regressionExample.predicted, '-s', ...
        'LineWidth', 1.1, 'MarkerSize', 4);
    hold off;
    xlim([1, 36]);
    xticks(1:5:36);
    xlabel('Division-1 part number');
    ylabel('LCModel concentration');
    title(sprintf('%s %s: observed vs PC1-PC3 prediction (adjusted R^2 = %.3f)', ...
        regressionExample.patientID, regressionExample.metabolite, ...
        regressionExample.adjustedR2), 'Interpreter', 'none');
    legend({'Observed', 'Predicted'}, 'Location', 'best');
    grid on;
    box on;
end


function [patientTable, summaryTable, globalSummary, example] = ...
    BuildPairResidualizationDiagnostics( ...
    pairwiseView, regressionSignals, examplePatientID, exampleMetabolites)

    sourceTable = pairwiseView.patientEligibilityTable;
    nRows = height(sourceTable);
    patientID = string(sourceTable.patientID);
    metaboliteA = string(sourceTable.metaboliteA);
    metaboliteB = string(sourceTable.metaboliteB);
    nValidPartsOriginal = double(sourceTable.nValidParts);
    nValidPartsCommon = zeros(nRows, 1);
    rEmpOriginal = nan(nRows, 1);
    rEmpResidual = nan(nRows, 1);
    rBaselineAssociated = nan(nRows, 1);
    rDeGraaf = nan(nRows, 1);
    deltaZOriginal = nan(nRows, 1);
    deltaZResidual = nan(nRows, 1);
    absDeltaZOriginal = nan(nRows, 1);
    absDeltaZResidual = nan(nRows, 1);
    improvementAbsDeltaZ = nan(nRows, 1);
    status = repmat("NOT_RUN", nRows, 1);
    baselineAssociatedStatus = repmat("NOT_RUN", nRows, 1);
    reason = strings(nRows, 1);

    signalPatientIDs = string({regressionSignals.patientID}).';
    for rowIndex = 1:nRows
        if ~sourceTable.eligible(rowIndex)
            status(rowIndex) = "ORIGINAL_INELIGIBLE";
            reason(rowIndex) = string(sourceTable.reason(rowIndex));
            continue;
        end

        signalIndex = find(signalPatientIDs == patientID(rowIndex), 1, 'first');
        if isempty(signalIndex)
            status(rowIndex) = "MISSING_PATIENT_SIGNALS";
            reason(rowIndex) = "Patient regression signals are unavailable.";
            continue;
        end
        signals = regressionSignals(signalIndex);
        metaboliteAIndex = find( ...
            signals.metaboliteList == metaboliteA(rowIndex), 1, 'first');
        metaboliteBIndex = find( ...
            signals.metaboliteList == metaboliteB(rowIndex), 1, 'first');
        if isempty(metaboliteAIndex) || isempty(metaboliteBIndex)
            status(rowIndex) = "MISSING_METABOLITE_SIGNALS";
            reason(rowIndex) = "One or both metabolite signals are unavailable.";
            continue;
        end
        if signals.regressionStatus(metaboliteAIndex) ~= "OK" || ...
                signals.regressionStatus(metaboliteBIndex) ~= "OK"
            status(rowIndex) = "REGRESSION_UNAVAILABLE";
            reason(rowIndex) = "One or both metabolite regressions are unavailable.";
            continue;
        end

        originalA = signals.concentration(:, metaboliteAIndex);
        originalB = signals.concentration(:, metaboliteBIndex);
        residualA = signals.residual(:, metaboliteAIndex);
        residualB = signals.residual(:, metaboliteBIndex);
        associatedA = signals.baselineAssociatedVariation(:, metaboliteAIndex);
        associatedB = signals.baselineAssociatedVariation(:, metaboliteBIndex);

        commonRows = isfinite(originalA) & isfinite(originalB) & ...
            isfinite(residualA) & isfinite(residualB);
        nValidPartsCommon(rowIndex) = sum(commonRows);
        if nValidPartsCommon(rowIndex) < pairwiseView.minValidParts
            status(rowIndex) = "INSUFFICIENT_COMMON_PARTS";
            reason(rowIndex) = sprintf( ...
                'Only %d common parts; minimum is %d.', ...
                nValidPartsCommon(rowIndex), pairwiseView.minValidParts);
            continue;
        end
        if nValidPartsCommon(rowIndex) ~= nValidPartsOriginal(rowIndex)
            error(['Patient %s pair %s/%s has %d original valid parts but ', ...
                '%d aligned residual parts.'], ...
                patientID(rowIndex), metaboliteA(rowIndex), metaboliteB(rowIndex), ...
                nValidPartsOriginal(rowIndex), nValidPartsCommon(rowIndex));
        end

        originalCheck = PearsonCorrelation( ...
            originalA(commonRows), originalB(commonRows));
        if ~isfinite(originalCheck) || ...
                abs(originalCheck - sourceTable.rEmpirical(rowIndex)) > 1e-10
            error(['Patient %s pair %s/%s does not reproduce the stored ', ...
                'original empirical correlation.'], ...
                patientID(rowIndex), metaboliteA(rowIndex), metaboliteB(rowIndex));
        end

        residualCorrelation = PearsonCorrelation( ...
            residualA(commonRows), residualB(commonRows));
        if ~isfinite(residualCorrelation)
            status(rowIndex) = "NONFINITE_RESIDUAL_CORRELATION";
            reason(rowIndex) = "Residual correlation is nonfinite.";
            continue;
        end

        associatedRows = commonRows & ...
            isfinite(associatedA) & isfinite(associatedB);
        associatedCorrelation = PearsonCorrelation( ...
            associatedA(associatedRows), associatedB(associatedRows));
        if isfinite(associatedCorrelation)
            baselineAssociatedStatus(rowIndex) = "OK";
        else
            baselineAssociatedStatus(rowIndex) = "ZERO_OR_NONFINITE_VARIANCE";
        end

        rEmpOriginal(rowIndex) = sourceTable.rEmpirical(rowIndex);
        rEmpResidual(rowIndex) = residualCorrelation;
        rBaselineAssociated(rowIndex) = associatedCorrelation;
        rDeGraaf(rowIndex) = sourceTable.rModel(rowIndex);
        zOriginal = atanh(ClampCorrelation(rEmpOriginal(rowIndex)));
        zResidual = atanh(ClampCorrelation(rEmpResidual(rowIndex)));
        zModel = atanh(ClampCorrelation(rDeGraaf(rowIndex)));
        deltaZOriginal(rowIndex) = zOriginal - zModel;
        deltaZResidual(rowIndex) = zResidual - zModel;
        absDeltaZOriginal(rowIndex) = abs(deltaZOriginal(rowIndex));
        absDeltaZResidual(rowIndex) = abs(deltaZResidual(rowIndex));
        improvementAbsDeltaZ(rowIndex) = ...
            absDeltaZOriginal(rowIndex) - absDeltaZResidual(rowIndex);
        status(rowIndex) = "OK";
        reason(rowIndex) = "";
    end

    patientTable = table( ...
        patientID, metaboliteA, metaboliteB, ...
        nValidPartsOriginal, nValidPartsCommon, ...
        rEmpOriginal, rEmpResidual, rBaselineAssociated, rDeGraaf, ...
        deltaZOriginal, deltaZResidual, ...
        absDeltaZOriginal, absDeltaZResidual, improvementAbsDeltaZ, ...
        status, baselineAssociatedStatus, reason);

    summaryTable = SummarizePairResidualization( ...
        patientTable, pairwiseView.pairTable, ...
        pairwiseView.minPatientsForGroupTest);
    validPairs = summaryTable.status == "OK";
    originalDistribution = summaryTable.meanDeltaZOriginal(validPairs);
    residualDistribution = summaryTable.meanDeltaZResidual(validPairs);
    if isempty(originalDistribution)
        error('No common metabolite pairs were available for histogram comparison.');
    end

    globalSummary = struct();
    globalSummary.nPairs = numel(originalDistribution);
    globalSummary.minPatientsRequired = ...
        pairwiseView.minPatientsForGroupTest;
    globalSummary.meanDeltaZOriginal = mean(originalDistribution);
    globalSummary.medianDeltaZOriginal = median(originalDistribution);
    globalSummary.sdDeltaZOriginal = std(originalDistribution, 0);
    globalSummary.meanDeltaZResidual = mean(residualDistribution);
    globalSummary.medianDeltaZResidual = median(residualDistribution);
    globalSummary.sdDeltaZResidual = std(residualDistribution, 0);
    globalSummary.absMeanShiftImprovement = ...
        abs(globalSummary.meanDeltaZOriginal) - ...
        abs(globalSummary.meanDeltaZResidual);
    globalSummary.fractionPairsCloserToZero = mean( ...
        abs(residualDistribution) < abs(originalDistribution));

    example = SelectPairResidualizationExample( ...
        patientTable, regressionSignals, examplePatientID, exampleMetabolites);
end


function summaryTable = SummarizePairResidualization( ...
    patientTable, sourcePairTable, minPatients)

    nPairs = height(sourcePairTable);
    metaboliteA = string(sourcePairTable.metaboliteA);
    metaboliteB = string(sourcePairTable.metaboliteB);
    nPatientsUsed = zeros(nPairs, 1);
    meanDeltaZOriginal = nan(nPairs, 1);
    medianDeltaZOriginal = nan(nPairs, 1);
    meanDeltaZResidual = nan(nPairs, 1);
    medianDeltaZResidual = nan(nPairs, 1);
    meanAbsDeltaZOriginal = nan(nPairs, 1);
    meanAbsDeltaZResidual = nan(nPairs, 1);
    meanImprovementAbsDeltaZ = nan(nPairs, 1);
    medianImprovementAbsDeltaZ = nan(nPairs, 1);
    meanREmpOriginal = nan(nPairs, 1);
    meanREmpResidual = nan(nPairs, 1);
    meanRBaselineAssociated = nan(nPairs, 1);
    meanRDeGraaf = nan(nPairs, 1);
    status = repmat("NOT_RUN", nPairs, 1);
    reason = strings(nPairs, 1);

    for pairIndex = 1:nPairs
        rows = patientTable.metaboliteA == metaboliteA(pairIndex) & ...
            patientTable.metaboliteB == metaboliteB(pairIndex) & ...
            patientTable.status == "OK";
        nPatientsUsed(pairIndex) = sum(rows);
        if nPatientsUsed(pairIndex) < minPatients
            status(pairIndex) = "INSUFFICIENT_PATIENTS";
            reason(pairIndex) = sprintf( ...
                'Only %d valid patients; minimum is %d.', ...
                nPatientsUsed(pairIndex), minPatients);
            continue;
        end

        original = patientTable.deltaZOriginal(rows);
        residual = patientTable.deltaZResidual(rows);
        improvement = patientTable.improvementAbsDeltaZ(rows);
        meanDeltaZOriginal(pairIndex) = mean(original);
        medianDeltaZOriginal(pairIndex) = median(original);
        meanDeltaZResidual(pairIndex) = mean(residual);
        medianDeltaZResidual(pairIndex) = median(residual);
        meanAbsDeltaZOriginal(pairIndex) = mean(abs(original));
        meanAbsDeltaZResidual(pairIndex) = mean(abs(residual));
        meanImprovementAbsDeltaZ(pairIndex) = mean(improvement);
        medianImprovementAbsDeltaZ(pairIndex) = median(improvement);
        meanREmpOriginal(pairIndex) = FisherMeanCorrelation( ...
            patientTable.rEmpOriginal(rows));
        meanREmpResidual(pairIndex) = FisherMeanCorrelation( ...
            patientTable.rEmpResidual(rows));
        meanRBaselineAssociated(pairIndex) = FisherMeanCorrelation( ...
            patientTable.rBaselineAssociated(rows));
        meanRDeGraaf(pairIndex) = FisherMeanCorrelation( ...
            patientTable.rDeGraaf(rows));
        status(pairIndex) = "OK";
        reason(pairIndex) = "";
    end

    summaryTable = table( ...
        metaboliteA, metaboliteB, nPatientsUsed, ...
        meanDeltaZOriginal, medianDeltaZOriginal, ...
        meanDeltaZResidual, medianDeltaZResidual, ...
        meanAbsDeltaZOriginal, meanAbsDeltaZResidual, ...
        meanImprovementAbsDeltaZ, medianImprovementAbsDeltaZ, ...
        meanREmpOriginal, meanREmpResidual, ...
        meanRBaselineAssociated, meanRDeGraaf, status, reason);
end


function example = SelectPairResidualizationExample( ...
    patientTable, regressionSignals, requestedPatientID, requestedMetabolites)

    rows = patientTable.status == "OK" & ...
        patientTable.patientID == requestedPatientID & ...
        ((patientTable.metaboliteA == requestedMetabolites(1) & ...
        patientTable.metaboliteB == requestedMetabolites(2)) | ...
        (patientTable.metaboliteA == requestedMetabolites(2) & ...
        patientTable.metaboliteB == requestedMetabolites(1)));
    rowIndex = find(rows, 1, 'first');
    if isempty(rowIndex)
        rows = patientTable.status == "OK" & ...
            patientTable.patientID == requestedPatientID;
        rowIndex = find(rows, 1, 'first');
    end
    if isempty(rowIndex)
        rowIndex = find(patientTable.status == "OK", 1, 'first');
    end
    if isempty(rowIndex)
        error('No valid patient-pair row is available for the example figure.');
    end

    patientID = patientTable.patientID(rowIndex);
    metaboliteA = patientTable.metaboliteA(rowIndex);
    metaboliteB = patientTable.metaboliteB(rowIndex);
    signalPatientIDs = string({regressionSignals.patientID}).';
    signalIndex = find(signalPatientIDs == patientID, 1, 'first');
    signals = regressionSignals(signalIndex);
    metaboliteAIndex = find(signals.metaboliteList == metaboliteA, 1, 'first');
    metaboliteBIndex = find(signals.metaboliteList == metaboliteB, 1, 'first');

    example = struct();
    example.patientID = patientID;
    example.metaboliteA = metaboliteA;
    example.metaboliteB = metaboliteB;
    example.partNumbers = signals.partNumbers;
    example.originalA = signals.concentration(:, metaboliteAIndex);
    example.predictedA = signals.predicted(:, metaboliteAIndex);
    example.residualA = signals.residual(:, metaboliteAIndex);
    example.originalB = signals.concentration(:, metaboliteBIndex);
    example.predictedB = signals.predicted(:, metaboliteBIndex);
    example.residualB = signals.residual(:, metaboliteBIndex);
    example.originalA(~isfinite(example.residualA)) = NaN;
    example.predictedA(~isfinite(example.residualA)) = NaN;
    example.originalB(~isfinite(example.residualB)) = NaN;
    example.predictedB(~isfinite(example.residualB)) = NaN;
end


function PrintPairResidualizationSummary(summaryTable, globalSummary)
    fprintf('\nBaseline-PC residual correlation diagnostic\n');
    fprintf('Common metabolite pairs: %d\n', globalSummary.nPairs);
    fprintf('Minimum patients per pair: %d\n', ...
        globalSummary.minPatientsRequired);
    fprintf('Original pair-level deltaZ: mean %.4f, median %.4f, SD %.4f\n', ...
        globalSummary.meanDeltaZOriginal, ...
        globalSummary.medianDeltaZOriginal, ...
        globalSummary.sdDeltaZOriginal);
    fprintf('Residual pair-level deltaZ: mean %.4f, median %.4f, SD %.4f\n', ...
        globalSummary.meanDeltaZResidual, ...
        globalSummary.medianDeltaZResidual, ...
        globalSummary.sdDeltaZResidual);
    fprintf('Absolute mean-shift improvement: %.4f\n', ...
        globalSummary.absMeanShiftImprovement);
    fprintf('Fraction of pairs moving closer to zero: %.4f\n', ...
        globalSummary.fractionPairsCloserToZero);

    validRows = find(summaryTable.status == "OK");
    [~, order] = sort(summaryTable.meanImprovementAbsDeltaZ(validRows), ...
        'descend');
    orderedRows = validRows(order);
    nShow = min(3, numel(orderedRows));
    fprintf('Largest mean absolute-deltaZ improvements:\n');
    for index = 1:nShow
        row = orderedRows(index);
        fprintf('  %s / %s: %.4f\n', ...
            summaryTable.metaboliteA(row), summaryTable.metaboliteB(row), ...
            summaryTable.meanImprovementAbsDeltaZ(row));
    end
    fprintf('Largest mean absolute-deltaZ worsenings:\n');
    for index = 1:nShow
        row = orderedRows(end - index + 1);
        fprintf('  %s / %s: %.4f\n', ...
            summaryTable.metaboliteA(row), summaryTable.metaboliteB(row), ...
            summaryTable.meanImprovementAbsDeltaZ(row));
    end
end


function PlotPairResidualizationFigures(summaryTable, example)
    validRows = summaryTable.status == "OK";
    original = summaryTable.meanDeltaZOriginal(validRows);
    residual = summaryTable.meanDeltaZResidual(validRows);
    maximumAbsolute = max(abs([original; residual]));
    if ~isfinite(maximumAbsolute) || maximumAbsolute <= 0
        maximumAbsolute = 1;
    end
    axisLimit = 1.05 * maximumAbsolute;
    binEdges = linspace(-axisLimit, axisLimit, 21);

    figure('Name', 'Baseline-adjusted empirical-minus-DeGraaf deltaZ', ...
        'NumberTitle', 'off');
    histogramLayout = tiledlayout(2, 1, ...
        'TileSpacing', 'compact', 'Padding', 'compact');
    nexttile;
    histogram(original, binEdges);
    xline(0, '--k');
    xlim([-axisLimit, axisLimit]);
    ylabel('Metabolite pairs');
    title('Original empirical minus DeGraaf Fisher-z difference');
    grid on;
    nexttile;
    histogram(residual, binEdges);
    xline(0, '--k');
    xlim([-axisLimit, axisLimit]);
    xlabel('Mean deltaZ across patients');
    ylabel('Metabolite pairs');
    title('After removing baseline-PC-associated concentration variation');
    grid on;
    title(histogramLayout, ...
        'Common-pair empirical-minus-DeGraaf deltaZ comparison');

    figure('Name', sprintf('%s %s-%s baseline residualization example', ...
        example.patientID, example.metaboliteA, example.metaboliteB), ...
        'NumberTitle', 'off');
    exampleLayout = tiledlayout(2, 1, ...
        'TileSpacing', 'compact', 'Padding', 'compact');
    nexttile;
    plot(example.partNumbers, example.originalA, '-o', 'LineWidth', 1.0);
    hold on;
    plot(example.partNumbers, example.predictedA, '-s', 'LineWidth', 1.0);
    plot(example.partNumbers, example.residualA, '-^', 'LineWidth', 1.0);
    hold off;
    xlim([1, 36]);
    ylabel(example.metaboliteA, 'Interpreter', 'none');
    title(sprintf('%s: original, PC prediction, and residual', ...
        example.metaboliteA), 'Interpreter', 'none');
    legend({'Original', 'PC prediction', 'Residual'}, 'Location', 'best');
    grid on;
    nexttile;
    plot(example.partNumbers, example.originalB, '-o', 'LineWidth', 1.0);
    hold on;
    plot(example.partNumbers, example.predictedB, '-s', 'LineWidth', 1.0);
    plot(example.partNumbers, example.residualB, '-^', 'LineWidth', 1.0);
    hold off;
    xlim([1, 36]);
    xlabel('Division-1 part number');
    ylabel(example.metaboliteB, 'Interpreter', 'none');
    title(sprintf('%s: original, PC prediction, and residual', ...
        example.metaboliteB), 'Interpreter', 'none');
    legend({'Original', 'PC prediction', 'Residual'}, 'Location', 'best');
    grid on;
    title(exampleLayout, sprintf('%s: %s / %s', ...
        example.patientID, example.metaboliteA, example.metaboliteB), ...
        'Interpreter', 'none');
end


function r = PearsonCorrelation(x, y)
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
    denominator = sqrt(sum(x.^2) .* sum(y.^2));
    if ~isfinite(denominator) || denominator <= 0
        r = NaN;
    else
        r = sum(x .* y) ./ denominator;
    end
end


function r = ClampCorrelation(r)
    r = max(min(r, 0.999999), -0.999999);
end


function r = FisherMeanCorrelation(values)
    values = double(values(:));
    values = values(isfinite(values));
    if isempty(values)
        r = NaN;
    else
        r = tanh(mean(atanh(ClampCorrelation(values))));
    end
end


function PrintCohortPCASummary(cohortTable)
    fprintf('\nCohort fitted-baseline PCA diagnostic\n');
    fprintf('Patients processed: %d\n', height(cohortTable));
    PrintMetricSummary('PC1 explained variance (%)', ...
        cohortTable.PC1ExplainedPct);
    PrintMetricSummary('Cumulative PC1-2 (%)', ...
        cohortTable.PC1To2CumulativePct);
    PrintMetricSummary('Cumulative PC1-3 (%)', ...
        cohortTable.PC1To3CumulativePct);
    PrintMetricSummary('Maximum PC1 energy fraction', ...
        cohortTable.maxPC1EnergyFraction);

    fprintf('Patients with cumulative PC1-2 >= 80%%: %d\n', ...
        sum(cohortTable.PC1To2CumulativePct >= 80));
    fprintf('Patients with cumulative PC1-2 >= 90%%: %d\n', ...
        sum(cohortTable.PC1To2CumulativePct >= 90));
    fprintf('Patients with cumulative PC1-3 >= 90%%: %d\n', ...
        sum(cohortTable.PC1To3CumulativePct >= 90));
    fprintf('Patients with cumulative PC1-3 >= 95%%: %d\n', ...
        sum(cohortTable.PC1To3CumulativePct >= 95));

    [maximumFraction, maximumIndex] = max( ...
        cohortTable.maxPC1EnergyFraction);
    fprintf('Largest maximum PC1 energy fraction: %.6f\n', maximumFraction);
    fprintf('  patient: %s\n', cohortTable.patientID(maximumIndex));
    fprintf('  part: %d\n', ...
        cohortTable.partWithLargestAbsPC1Score(maximumIndex));
end


function PrintMetricSummary(label, values)
    fprintf('%s: mean %.2f, median %.2f, min %.2f, max %.2f\n', ...
        label, mean(values), median(values), min(values), max(values));
end


function PlotCohortPCAFigures(cohortTable)
    cumulativeVariance = [ ...
        cohortTable.PC1CumulativePct, ...
        cohortTable.PC1To2CumulativePct, ...
        cohortTable.PC1To3CumulativePct];

    figure('Name', 'Cohort baseline PCA cumulative explained variance', ...
        'NumberTitle', 'off');
    if height(cohortTable) == 1
        bar(1:3, cumulativeVariance);
        xticks(1:3);
        xticklabels({'PC1', 'PC1-2', 'PC1-3'});
    else
        boxplot(cumulativeVariance, ...
            'Labels', {'PC1', 'PC1-2', 'PC1-3'}, ...
            'Symbol', '');
        hold on;
        patientOffsets = linspace(-0.12, 0.12, height(cohortTable)).';
        for cumulativeIndex = 1:3
            scatter(cumulativeIndex + patientOffsets, ...
                cumulativeVariance(:, cumulativeIndex), 14, ...
                [0.25, 0.25, 0.25], 'filled', ...
                'MarkerFaceAlpha', 0.55);
        end
        hold off;
    end
    ylim([0, 100]);
    ylabel('Cumulative explained variance (%)');
    title('Patient-specific fitted-baseline PCA variance');
    grid on;
    box on;

    figure('Name', 'Cohort baseline PCA PC1 dominance QC', ...
        'NumberTitle', 'off');
    patientIndex = (1:height(cohortTable)).';
    plot(patientIndex, cohortTable.maxPC1EnergyFraction, '-o', ...
        'LineWidth', 1.0, 'MarkerSize', 4);
    if height(cohortTable) == 1
        xlim([0.5, 1.5]);
    else
        xlim([1, height(cohortTable)]);
    end
    xticks(patientIndex);
    xticklabels(cohortTable.patientID);
    xtickangle(90);
    set(gca, 'TickLabelInterpreter', 'none');
    xlabel('Patient');
    ylabel('Maximum PC1 score-energy fraction');
    title('PC1 dominance by a single Division-1 part');
    grid on;
    box on;
end


function [baselineMatrix, ppmAxis] = ExtractDivision1Baselines( ...
    patientID, coordFiles, fitData)

    coordFiles = string(coordFiles(:));
    fitData = fitData(:);

    if numel(coordFiles) ~= numel(fitData)
        error(['Patient %s has %d coord filenames but %d fitData entries; ', ...
            'their source correspondence cannot be verified.'], ...
            patientID, numel(coordFiles), numel(fitData));
    end

    division1Indices = zeros(36, 1);

    for fileIndex = 1:numel(coordFiles)
        [~, filename, extension] = fileparts(coordFiles(fileIndex));
        tokens = regexp(filename + extension, ...
            '(^|_)Division_1_part_(\d+)\.basis\.coord$', ...
            'tokens', 'once');

        if isempty(tokens)
            continue;
        end

        partNumber = str2double(tokens{2});
        if partNumber < 1 || partNumber > 36 || partNumber ~= fix(partNumber)
            error('Patient %s has unexpected Division-1 part number %g in %s.', ...
                patientID, partNumber, coordFiles(fileIndex));
        end
        if division1Indices(partNumber) ~= 0
            error('Patient %s has duplicate Division-1 part %d.', ...
                patientID, partNumber);
        end

        division1Indices(partNumber) = fileIndex;
    end

    missingParts = find(division1Indices == 0);
    if ~isempty(missingParts)
        error('Patient %s is missing Division-1 part(s): %s.', ...
            patientID, strjoin(string(missingParts), ', '));
    end

    baselineMatrix = nan(36, numel(fitData(division1Indices(1)).baseline));
    ppmAxis = [];

    for partNumber = 1:36
        fitEntry = fitData(division1Indices(partNumber));

        if ~isfield(fitEntry, 'baseline') || isempty(fitEntry.baseline) || ...
                ~isnumeric(fitEntry.baseline) || ~isvector(fitEntry.baseline) || ...
                any(~isfinite(fitEntry.baseline(:)))
            error(['Patient %s Division-1 part %d has a baseline that is ', ...
                'missing, empty, nonnumeric, nonvector, or nonfinite.'], ...
                patientID, partNumber);
        end
        if ~isfield(fitEntry, 'axis') || isempty(fitEntry.axis) || ...
                ~isnumeric(fitEntry.axis) || ~isvector(fitEntry.axis) || ...
                any(~isfinite(fitEntry.axis(:)))
            error(['Patient %s Division-1 part %d has a ppm grid that is ', ...
                'missing, empty, nonnumeric, nonvector, or nonfinite.'], ...
                patientID, partNumber);
        end

        baseline = reshape(fitEntry.baseline, 1, []);
        currentPPM = reshape(fitEntry.axis, 1, []);

        if numel(baseline) ~= numel(currentPPM)
            error(['Patient %s Division-1 part %d baseline length (%d) ', ...
                'does not match ppm-grid length (%d).'], ...
                patientID, partNumber, numel(baseline), numel(currentPPM));
        end

        if partNumber == 1
            ppmAxis = currentPPM;
        elseif ~isequal(currentPPM, ppmAxis)
            error(['Patient %s Division-1 part %d does not use the same ', ...
                'ppm grid and ordering as part 1.'], patientID, partNumber);
        end

        baselineMatrix(partNumber, :) = baseline;
    end
end
