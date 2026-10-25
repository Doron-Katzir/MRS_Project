%% TwoGaussianT2CovarianceMonteCarlo
% Monte Carlo study of covariance-corrected T2 estimation for two
% overlapping Gaussian spectral basis functions.
%
% At each echo time, the K-by-1 spectrum is
%
%   s_i = B(d)*a_i + epsilon_i,
%   a_i = [A1*exp(-TE_i/T2_1); A2*exp(-TE_i/T2_2)],
%
% with independent epsilon_i ~ N(0, sigmaNoise^2*I_K). The experiment has
% no random concentration effects and no random T2 parameters.
%
% The script asks two distinct questions:
%   1. Does subtracting the known OLS amplitude-noise covariance improve
%      estimation of the covariance of the noise-free decay trajectories?
%   2. Does any covariance improvement translate to improved T2 estimates?

clear;
clc;
close all;

%% Configuration (all simulation and numerical settings are here)

cfg.K = 160;
cfg.x = linspace(0, 10, cfg.K).';
cfg.gaussianWidth = 0.55;
cfg.separationValues = linspace(0.2, 3, 15).';
cfg.trueAmplitude = [1, 1];
cfg.echoTimesMs = (10:20:190).';
cfg.sigmaNoise = 0.05;
cfg.nMonteCarlo = 500;
cfg.rngSeed = 1;

cfg.scenarioNames = ["Equal T2 (100, 100 ms)", ...
    "Different T2 (100, 300 ms)"];
cfg.trueT2Ms = [100, 100; 100, 300];

% Conventional exponential-NLS settings. These fixed initial values are
% not equal to either simulated T2 parameter set.
cfg.nlsInitialAmplitude = [0.8, 0.8];
cfg.nlsInitialT2Ms = [180, 180];
cfg.nlsAmplitudeLowerBound = [0, 0];
cfg.nlsAmplitudeUpperBound = [5, 5];
cfg.t2LowerBoundMs = 10;
cfg.t2UpperBoundMs = 2000;
cfg.boundHitRelativeTolerance = 1e-6;

cfg.optimizer.functionTolerance = 1e-10;
cfg.optimizer.stepTolerance = 1e-10;
cfg.optimizer.optimalityTolerance = 1e-8;
cfg.optimizer.maxIterations = 500;
cfg.optimizer.maxFunctionEvaluations = 4000;
cfg.optimizer.display = "off";
cfg.rootToleranceMs = 1e-9;

% Validation thresholds express numerical precision or sampling
% uncertainty; they do not encode expected scientific outcomes.
cfg.validation.projectionTolerance = 1e-10;
cfg.validation.noiseCovarianceTolerance = 1e-12;
cfg.validation.ratioT2ToleranceMs = 1e-6;
cfg.validation.standardErrorMultiplier = 6;
cfg.validation.minimumMonteCarloForSamplingCheck = 30;

cfg.representativeSeparation = 0.4;
cfg.makeFigures = true;
cfg.saveOutputs = true;
cfg.figureVisible = "off";
cfg.outputDirectory = fullfile(fileparts(mfilename("fullpath")), ...
    "T2CovarianceMonteCarloResults");

% Set MRS_T2_COV_SMOKE_TEST=1 only for a fast automated code check. The
% scientific defaults above remain unchanged for ordinary execution.
if strcmp(getenv("MRS_T2_COV_SMOKE_TEST"), "1")
    cfg.nMonteCarlo = 3;
    cfg.separationValues = cfg.separationValues(1:2);
    cfg.makeFigures = false;
    cfg.saveOutputs = false;
end

validateConfiguration(cfg);
rng(cfg.rngSeed, "twister");

t2MethodNames = ["Uncorrected covariance", ...
    "Corrected covariance", "Conventional exponential NLS"];
covarianceMethodNames = ["Uncorrected OLS covariance", ...
    "Corrected OLS covariance", "Model-based NLS covariance"];
nT2Methods = numel(t2MethodNames);
nCovarianceMethods = numel(covarianceMethodNames);
nScenarios = size(cfg.trueT2Ms, 1);
nSeparations = numel(cfg.separationValues);
nEchoes = numel(cfg.echoTimesMs);

%% Precompute bases, OLS projections, and analytical noise covariances

basisCache = repmat(struct( ...
    "B", [], "projection", [], "noiseCovariance", [], ...
    "basisCorrelation", NaN, "gramConditionNumber", NaN), ...
    nSeparations, 1);

for dIndex = 1:nSeparations
    d = cfg.separationValues(dIndex);
    centers = [5 - d/2, 5 + d/2];
    b1 = exp(-0.5*((cfg.x - centers(1))/cfg.gaussianWidth).^2);
    b2 = exp(-0.5*((cfg.x - centers(2))/cfg.gaussianWidth).^2);
    b1 = b1/max(b1);
    b2 = b2/max(b2);
    B = [b1, b2];
    gram = B.'*B;

    % Both the OLS projection and inverse Gram factor are obtained by
    % linear solves; no explicit matrix inversion is used.
    basisCache(dIndex).B = B;
    basisCache(dIndex).projection = B\eye(cfg.K);
    basisCache(dIndex).noiseCovariance = ...
        cfg.sigmaNoise^2*(gram\eye(2));
    basisCache(dIndex).basisCorrelation = ...
        (b1.'*b2)/sqrt((b1.'*b1)*(b2.'*b2));
    basisCache(dIndex).gramConditionNumber = cond(gram);
end

%% Deterministic true trajectories and covariance references

trueTrajectories = nan(2, nEchoes, nScenarios);
trueCovariance = nan(2, 2, nScenarios);
for scenarioIndex = 1:nScenarios
    trueTrajectories(:, :, scenarioIndex) = ...
        cfg.trueAmplitude(:).*exp( ...
        -cfg.echoTimesMs.'./cfg.trueT2Ms(scenarioIndex, :).');
    trueCovariance(:, :, scenarioIndex) = sampleCovariance( ...
        trueTrajectories(:, :, scenarioIndex));
end

%% Mathematical validation checks available before Monte Carlo sampling

preValidation = runPreSimulationValidation( ...
    cfg, basisCache, trueTrajectories);
assert(all(preValidation.projectionPassed, "all"), ...
    "Noise-free OLS projection validation failed.");
assert(all(preValidation.noiseCovariancePassed), ...
    "Analytical OLS noise-covariance validation failed.");
assert(all(preValidation.ratioRecoveryPassed, "all"), ...
    "Noise-free T2-ratio inversion validation failed.");
assert(preValidation.ratioMonotonicOnBounds, ...
    "The theoretical T2 ratio is not monotone on the configured bounds.");

%% Allocate raw Monte Carlo results

raw = struct();
raw.t2MethodNames = t2MethodNames;
raw.covarianceMethodNames = covarianceMethodNames;
raw.scenarioNames = cfg.scenarioNames;
raw.separationValues = cfg.separationValues;
raw.basisCorrelation = reshape([basisCache.basisCorrelation], [], 1);
raw.gramConditionNumber = reshape( ...
    [basisCache.gramConditionNumber], [], 1);
raw.trueTrajectories = trueTrajectories;
raw.trueCovariance = trueCovariance;

raw.olsAmplitudeHat = nan( ...
    2, nEchoes, nScenarios, nSeparations, cfg.nMonteCarlo);
raw.olsMeanAmplitude = nan( ...
    2, nScenarios, nSeparations, cfg.nMonteCarlo);
raw.nlsFittedTrajectory = nan( ...
    2, nEchoes, nScenarios, nSeparations, cfg.nMonteCarlo);

% Covariance layout: 2-by-2-by-scenario-by-separation-by-MC-by-method.
raw.covarianceEstimate = nan(2, 2, nScenarios, nSeparations, ...
    cfg.nMonteCarlo, nCovarianceMethods);
raw.frobeniusError = nan( ...
    nScenarios, nSeparations, cfg.nMonteCarlo, nCovarianceMethods);
raw.correctedHasNegativeDiagonal = false( ...
    nScenarios, nSeparations, cfg.nMonteCarlo);
raw.correctedMinEigenvalue = nan( ...
    nScenarios, nSeparations, cfg.nMonteCarlo);

% T2 layout: scenario-by-separation-by-MC-by-method-by-metabolite.
t2ResultSize = [nScenarios, nSeparations, cfg.nMonteCarlo, ...
    nT2Methods, 2];
raw.t2HatMs = nan(t2ResultSize);
raw.t2EstimateValid = false(t2ResultSize);
raw.ratioHat = nan(t2ResultSize);
raw.t2StatusCode = zeros(t2ResultSize, "uint8");
raw.t2BoundHit = nan(t2ResultSize);

% Ratio status codes: 0 not evaluated, 1 valid, 2 negative ratio,
% 3 nonfinite ratio/result, 4 outside the theoretical range, and
% 5 bounded root-finding failure. NLS uses codes 1, 3, 5, and 6, where
% 6 denotes a nonpositive optimizer exit flag.
raw.nlsAmplitudeHat = nan( ...
    nScenarios, nSeparations, cfg.nMonteCarlo, 2);
raw.nlsCandidateAmplitude = nan( ...
    nScenarios, nSeparations, cfg.nMonteCarlo, 2);
raw.nlsCandidateT2Ms = nan( ...
    nScenarios, nSeparations, cfg.nMonteCarlo, 2);
raw.nlsOptimizerConverged = false( ...
    nScenarios, nSeparations, cfg.nMonteCarlo, 2);
raw.nlsExitFlag = nan( ...
    nScenarios, nSeparations, cfg.nMonteCarlo, 2);
raw.nlsIterations = nan( ...
    nScenarios, nSeparations, cfg.nMonteCarlo, 2);
raw.nlsFunctionCount = nan( ...
    nScenarios, nSeparations, cfg.nMonteCarlo, 2);
raw.nlsAmplitudeBoundHit = false( ...
    nScenarios, nSeparations, cfg.nMonteCarlo, 2);
raw.nlsOptimizerName = strings( ...
    nScenarios, nSeparations, cfg.nMonteCarlo, 2);

%% Monte Carlo experiment

for scenarioIndex = 1:nScenarios
    noiseFreeAmplitude = trueTrajectories(:, :, scenarioIndex);
    referenceCovariance = trueCovariance(:, :, scenarioIndex);

    for dIndex = 1:nSeparations
        cache = basisCache(dIndex);
        cleanSpectra = cache.B*noiseFreeAmplitude;

        fprintf("Scenario %d/%d, d = %.3f (%d/%d), cond(B''B) = %.3g\n", ...
            scenarioIndex, nScenarios, cfg.separationValues(dIndex), ...
            dIndex, nSeparations, cache.gramConditionNumber);

        for mcIndex = 1:cfg.nMonteCarlo
            spectra = cleanSpectra + ...
                cfg.sigmaNoise*randn(cfg.K, nEchoes);
            amplitudeOLS = cache.projection*spectra;
            meanAmplitudeOLS = mean(amplitudeOLS, 2);
            covarianceOLS = sampleCovariance(amplitudeOLS);
            covarianceCorrected = covarianceOLS - cache.noiseCovariance;

            raw.olsAmplitudeHat(:, :, scenarioIndex, dIndex, mcIndex) = ...
                amplitudeOLS;
            raw.olsMeanAmplitude(:, scenarioIndex, dIndex, mcIndex) = ...
                meanAmplitudeOLS;
            raw.covarianceEstimate(:, :, scenarioIndex, dIndex, ...
                mcIndex, 1) = covarianceOLS;
            raw.covarianceEstimate(:, :, scenarioIndex, dIndex, ...
                mcIndex, 2) = covarianceCorrected;

            correctedEigenvalues = eig((covarianceCorrected + ...
                covarianceCorrected.')/2);
            raw.correctedHasNegativeDiagonal(scenarioIndex, dIndex, ...
                mcIndex) = any(diag(covarianceCorrected) < 0);
            raw.correctedMinEigenvalue(scenarioIndex, dIndex, mcIndex) = ...
                min(correctedEigenvalues);

            % Methods A and B use the identical bounded ratio-inversion
            % function; their only difference is the covariance diagonal.
            covariancePair = cat(3, covarianceOLS, covarianceCorrected);
            for methodIndex = 1:2
                for metaboliteIndex = 1:2
                    ratioHat = covariancePair(metaboliteIndex, ...
                        metaboliteIndex, methodIndex)/ ...
                        meanAmplitudeOLS(metaboliteIndex)^2;
                    ratioFit = invertT2Ratio(ratioHat, cfg.echoTimesMs, ...
                        cfg.t2LowerBoundMs, cfg.t2UpperBoundMs, ...
                        cfg.rootToleranceMs, ...
                        cfg.boundHitRelativeTolerance);
                    raw.ratioHat(scenarioIndex, dIndex, mcIndex, ...
                        methodIndex, metaboliteIndex) = ratioHat;
                    raw.t2StatusCode(scenarioIndex, dIndex, mcIndex, ...
                        methodIndex, metaboliteIndex) = ratioFit.statusCode;
                    raw.t2BoundHit(scenarioIndex, dIndex, mcIndex, ...
                        methodIndex, metaboliteIndex) = ...
                        double(ratioFit.boundHit);
                    if ratioFit.valid
                        raw.t2HatMs(scenarioIndex, dIndex, mcIndex, ...
                            methodIndex, metaboliteIndex) = ratioFit.t2Ms;
                        raw.t2EstimateValid(scenarioIndex, dIndex, mcIndex, ...
                            methodIndex, metaboliteIndex) = true;
                    end
                end
            end

            % Method C fits each OLS amplitude trajectory conventionally.
            nlsSuccess = false(1, 2);
            fittedTrajectories = nan(2, nEchoes);
            for metaboliteIndex = 1:2
                residualFunction = @(parameter) parameter(1)*exp( ...
                    -cfg.echoTimesMs/parameter(2)) - ...
                    amplitudeOLS(metaboliteIndex, :).';
                nlsFit = solveBoundedLeastSquares(residualFunction, ...
                    [cfg.nlsInitialAmplitude(metaboliteIndex), ...
                    cfg.nlsInitialT2Ms(metaboliteIndex)], ...
                    [cfg.nlsAmplitudeLowerBound(metaboliteIndex), ...
                    cfg.t2LowerBoundMs], ...
                    [cfg.nlsAmplitudeUpperBound(metaboliteIndex), ...
                    cfg.t2UpperBoundMs], cfg.optimizer, ...
                    cfg.boundHitRelativeTolerance);

                raw.nlsCandidateAmplitude(scenarioIndex, dIndex, mcIndex, ...
                    metaboliteIndex) = nlsFit.candidate(1);
                raw.nlsCandidateT2Ms(scenarioIndex, dIndex, mcIndex, ...
                    metaboliteIndex) = nlsFit.candidate(2);
                raw.nlsOptimizerConverged(scenarioIndex, dIndex, mcIndex, ...
                    metaboliteIndex) = nlsFit.success;
                raw.nlsExitFlag(scenarioIndex, dIndex, mcIndex, ...
                    metaboliteIndex) = nlsFit.exitFlag;
                raw.nlsIterations(scenarioIndex, dIndex, mcIndex, ...
                    metaboliteIndex) = nlsFit.iterations;
                raw.nlsFunctionCount(scenarioIndex, dIndex, mcIndex, ...
                    metaboliteIndex) = nlsFit.functionCount;
                raw.nlsAmplitudeBoundHit(scenarioIndex, dIndex, mcIndex, ...
                    metaboliteIndex) = nlsFit.boundHit(1);
                raw.nlsOptimizerName(scenarioIndex, dIndex, mcIndex, ...
                    metaboliteIndex) = nlsFit.optimizerName;
                raw.t2BoundHit(scenarioIndex, dIndex, mcIndex, 3, ...
                    metaboliteIndex) = double(nlsFit.boundHit(2));

                if nlsFit.success
                    nlsSuccess(metaboliteIndex) = true;
                    raw.nlsAmplitudeHat(scenarioIndex, dIndex, mcIndex, ...
                        metaboliteIndex) = nlsFit.candidate(1);
                    raw.t2HatMs(scenarioIndex, dIndex, mcIndex, 3, ...
                        metaboliteIndex) = nlsFit.candidate(2);
                    raw.t2EstimateValid(scenarioIndex, dIndex, mcIndex, 3, ...
                        metaboliteIndex) = true;
                    raw.t2StatusCode(scenarioIndex, dIndex, mcIndex, 3, ...
                        metaboliteIndex) = uint8(1);
                    fittedTrajectories(metaboliteIndex, :) = ...
                        nlsFit.candidate(1)*exp( ...
                        -cfg.echoTimesMs.'/nlsFit.candidate(2));
                else
                    raw.t2StatusCode(scenarioIndex, dIndex, mcIndex, 3, ...
                        metaboliteIndex) = nlsFit.statusCode;
                end
            end

            raw.nlsFittedTrajectory(:, :, scenarioIndex, dIndex, mcIndex) = ...
                fittedTrajectories;
            if all(nlsSuccess)
                covarianceNLS = sampleCovariance(fittedTrajectories);
                raw.covarianceEstimate(:, :, scenarioIndex, dIndex, ...
                    mcIndex, 3) = covarianceNLS;
            end

            for covarianceMethodIndex = 1:nCovarianceMethods
                covarianceEstimate = raw.covarianceEstimate(:, :, ...
                    scenarioIndex, dIndex, mcIndex, covarianceMethodIndex);
                if all(isfinite(covarianceEstimate), "all")
                    raw.frobeniusError(scenarioIndex, dIndex, mcIndex, ...
                        covarianceMethodIndex) = norm( ...
                        covarianceEstimate - referenceCovariance, "fro");
                end
            end
        end
    end
end

%% Summaries and empirical covariance-bias validation

t2SummaryTable = summarizeT2Results(raw, cfg);
[covarianceSummaryTable, covarianceStats] = ...
    summarizeCovarianceResults(raw, cfg);
validityTable = summarizeValidity(raw, cfg);
covarianceDiagnosticsTable = summarizeCorrectedCovarianceDiagnostics( ...
    raw, cfg);
samplingValidationTable = validateExpectedCovarianceBias( ...
    raw, cfg, basisCache);

samplingChecksEvaluated = cfg.nMonteCarlo >= ...
    cfg.validation.minimumMonteCarloForSamplingCheck;
if samplingChecksEvaluated
    covarianceBiasValidationPassed = ...
        all(samplingValidationTable.Passed);
    if ~covarianceBiasValidationPassed
        warning("TwoGaussianT2CovarianceMonteCarlo:SamplingValidation", ...
            "At least one covariance-mean check exceeded %.1f Monte " + ...
            "Carlo standard errors.", ...
            cfg.validation.standardErrorMultiplier);
    end
else
    covarianceBiasValidationPassed = NaN;
end

validationResults = struct();
validationResults.preSimulation = preValidation;
validationResults.samplingTable = samplingValidationTable;
validationResults.samplingChecksEvaluated = samplingChecksEvaluated;
validationResults.covarianceBiasValidationPassed = ...
    covarianceBiasValidationPassed;

disp("T2 metrics (bias/RMSE/median error are conditional on validity):");
disp(t2SummaryTable);
disp("Covariance metrics (all methods use a common valid repetition set):");
disp(covarianceSummaryTable);
disp("Validity and optimizer diagnostics:");
disp(validityTable);
disp("Corrected-covariance indefiniteness diagnostics:");
disp(covarianceDiagnosticsTable);
disp("Expected covariance-bias sampling checks:");
disp(samplingValidationTable);

%% Save data tables before figure generation

if cfg.saveOutputs
    if ~isfolder(cfg.outputDirectory)
        mkdir(cfg.outputDirectory);
    end
    writetable(t2SummaryTable, fullfile(cfg.outputDirectory, ...
        "TwoGaussianT2Covariance_T2Summary.csv"));
    writetable(covarianceSummaryTable, fullfile(cfg.outputDirectory, ...
        "TwoGaussianT2Covariance_CovarianceSummary.csv"));
    writetable(validityTable, fullfile(cfg.outputDirectory, ...
        "TwoGaussianT2Covariance_Validity.csv"));
    writetable(covarianceDiagnosticsTable, fullfile(cfg.outputDirectory, ...
        "TwoGaussianT2Covariance_Indefiniteness.csv"));
    writetable(samplingValidationTable, fullfile(cfg.outputDirectory, ...
        "TwoGaussianT2Covariance_Validation.csv"));
    save(fullfile(cfg.outputDirectory, ...
        "TwoGaussianT2Covariance_RawResults.mat"), ...
        "raw", "cfg", "basisCache", "t2SummaryTable", ...
        "covarianceSummaryTable", "covarianceStats", ...
        "validityTable", "covarianceDiagnosticsTable", ...
        "validationResults", "-v7.3");
end

%% Figure 1: conditional-on-valid T2 RMSE versus separation

if cfg.makeFigures
    colors = lines(nT2Methods);
    markers = ["o", "s", "^"];

    fig1 = figure("Color", "w", "Visible", cfg.figureVisible, ...
        "Name", "T2 RMSE versus peak separation");
    layout = tiledlayout(fig1, 2, 2, "TileSpacing", "compact", ...
        "Padding", "compact");
    for scenarioIndex = 1:nScenarios
        for metaboliteIndex = 1:2
            ax = nexttile(layout, (scenarioIndex - 1)*2 + metaboliteIndex);
            hold(ax, "on");
            for methodIndex = 1:nT2Methods
                rows = t2SummaryTable.ScenarioIndex == scenarioIndex & ...
                    t2SummaryTable.Method == t2MethodNames(methodIndex);
                y = t2SummaryTable.(sprintf( ...
                    "RMSEConditional_T2_%d_ms", metaboliteIndex))(rows);
                plot(ax, t2SummaryTable.Separation(rows), y, ...
                    "-" + markers(methodIndex), ...
                    "Color", colors(methodIndex, :), ...
                    "LineWidth", 1.6, "MarkerSize", 6, ...
                    "DisplayName", t2MethodNames(methodIndex));
            end
            xlabel(ax, "Gaussian center separation d");
            ylabel(ax, sprintf("T2_%d RMSE (ms)", metaboliteIndex));
            title(ax, sprintf("%s; T2_%d = %.0f ms", ...
                cfg.scenarioNames(scenarioIndex), metaboliteIndex, ...
                cfg.trueT2Ms(scenarioIndex, metaboliteIndex)), ...
                "Interpreter", "none");
            grid(ax, "on");
            stylePlotAxes(ax);
            if scenarioIndex == 1 && metaboliteIndex == 1
                lgd = legend(ax, "Location", "best");
                set(lgd, "Color", "w", "TextColor", "k", ...
                    "EdgeColor", [0.3, 0.3, 0.3]);
            end
        end
    end
    layoutTitle = title(layout, ...
        "T2 RMSE conditional on valid estimates");
    layoutTitle.Color = "k";
    saveFigureFiles(fig1, cfg, "Figure1_T2_RMSE");

    %% Figure 2: conditional-on-valid T2 bias versus separation

    fig2 = figure("Color", "w", "Visible", cfg.figureVisible, ...
        "Name", "T2 bias versus peak separation");
    layout = tiledlayout(fig2, 2, 2, "TileSpacing", "compact", ...
        "Padding", "compact");
    for scenarioIndex = 1:nScenarios
        for metaboliteIndex = 1:2
            ax = nexttile(layout, (scenarioIndex - 1)*2 + metaboliteIndex);
            hold(ax, "on");
            yline(ax, 0, "k--", "HandleVisibility", "off");
            for methodIndex = 1:nT2Methods
                rows = t2SummaryTable.ScenarioIndex == scenarioIndex & ...
                    t2SummaryTable.Method == t2MethodNames(methodIndex);
                y = t2SummaryTable.(sprintf( ...
                    "BiasConditional_T2_%d_ms", metaboliteIndex))(rows);
                plot(ax, t2SummaryTable.Separation(rows), y, ...
                    "-" + markers(methodIndex), ...
                    "Color", colors(methodIndex, :), ...
                    "LineWidth", 1.6, "MarkerSize", 6, ...
                    "DisplayName", t2MethodNames(methodIndex));
            end
            xlabel(ax, "Gaussian center separation d");
            ylabel(ax, sprintf("T2_%d bias (ms)", metaboliteIndex));
            title(ax, cfg.scenarioNames(scenarioIndex), ...
                "Interpreter", "none");
            grid(ax, "on");
            stylePlotAxes(ax);
            if scenarioIndex == 1 && metaboliteIndex == 1
                lgd = legend(ax, "Location", "best");
                set(lgd, "Color", "w", "TextColor", "k", ...
                    "EdgeColor", [0.3, 0.3, 0.3]);
            end
        end
    end
    layoutTitle = title(layout, "T2 bias conditional on valid estimates");
    layoutTitle.Color = "k";
    saveFigureFiles(fig2, cfg, "Figure2_T2_Bias");

    %% Figure 3: covariance root mean squared Frobenius error

    fig3 = figure("Color", "w", "Visible", cfg.figureVisible, ...
        "Name", "Covariance error versus peak separation");
    layout = tiledlayout(fig3, 1, 2, "TileSpacing", "compact", ...
        "Padding", "compact");
    for scenarioIndex = 1:nScenarios
        ax = nexttile(layout);
        hold(ax, "on");
        for methodIndex = 1:nCovarianceMethods
            rows = covarianceSummaryTable.ScenarioIndex == scenarioIndex & ...
                covarianceSummaryTable.Method == ...
                covarianceMethodNames(methodIndex);
            plot(ax, covarianceSummaryTable.Separation(rows), ...
                covarianceSummaryTable.RMSE_Frobenius(rows), ...
                "-" + markers(methodIndex), ...
                "Color", colors(methodIndex, :), ...
                "LineWidth", 1.6, "MarkerSize", 6, ...
                "DisplayName", covarianceMethodNames(methodIndex));
        end
        xlabel(ax, "Gaussian center separation d");
        ylabel(ax, "Covariance RMSE (Frobenius norm)");
        title(ax, cfg.scenarioNames(scenarioIndex), ...
            "Interpreter", "none");
        grid(ax, "on");
        stylePlotAxes(ax);
        if scenarioIndex == 1
            lgd = legend(ax, "Location", "best");
            set(lgd, "Color", "w", "TextColor", "k", ...
                "EdgeColor", [0.3, 0.3, 0.3]);
        end
    end
    layoutTitle = title(layout, "Covariance estimation error");
    layoutTitle.Color = "k";
    saveFigureFiles(fig3, cfg, "Figure3_Covariance_RMSE");

    %% Figure 4: covariance matrices at representative separation

    [~, representativeIndex] = min(abs( ...
        cfg.separationValues - cfg.representativeSeparation));
    representativeD = cfg.separationValues(representativeIndex);
    fig4 = figure("Color", "w", "Visible", cfg.figureVisible, ...
        "Name", "Covariance matrix comparison");
    layout = tiledlayout(fig4, 2, 4, "TileSpacing", "compact", ...
        "Padding", "compact");
    matrixTitles = ["True covariance", ...
        "Mean uncorrected", "Mean corrected", "Mean model NLS"];

    for scenarioIndex = 1:nScenarios
        matrixSet = nan(2, 2, 4);
        matrixSet(:, :, 1) = trueCovariance(:, :, scenarioIndex);
        matrixSet(:, :, 2:4) = covarianceStats.meanEstimate(:, :, ...
            scenarioIndex, representativeIndex, :);
        colorLimits = [min(matrixSet, [], "all"), ...
            max(matrixSet, [], "all")];
        if diff(colorLimits) <= eps(max(1, max(abs(colorLimits))))
            padding = max(1, abs(colorLimits(1)))*1e-6;
            colorLimits = colorLimits + [-padding, padding];
        end

        for matrixIndex = 1:4
            ax = nexttile(layout, (scenarioIndex - 1)*4 + matrixIndex);
            imagesc(ax, matrixSet(:, :, matrixIndex), colorLimits);
            axis(ax, "image");
            xticks(ax, 1:2);
            yticks(ax, 1:2);
            xticklabels(ax, ["Metabolite 1", "Metabolite 2"]);
            if matrixIndex == 1
                yticklabels(ax, ["Metabolite 1", "Metabolite 2"]);
            else
                yticklabels(ax, []);
            end
            xtickangle(ax, 25);
            title(ax, matrixTitles(matrixIndex), "Interpreter", "none");
            addMatrixValueLabels(ax, matrixSet(:, :, matrixIndex), ...
                colorLimits);
            if matrixIndex == 1
                ylabel(ax, cfg.scenarioNames(scenarioIndex), ...
                    "Interpreter", "none", "FontWeight", "bold");
            end
            stylePlotAxes(ax);
            if matrixIndex == 4
                cb = colorbar(ax);
                cb.Color = "k";
            end
        end
    end
    colormap(fig4, parula(256));
    layoutTitle = title(layout, sprintf( ...
        "Covariance comparison at d = %.3g", representativeD));
    layoutTitle.Color = "k";
    saveFigureFiles(fig4, cfg, "Figure4_Covariance_Matrices");
end

%% Local functions

function validateConfiguration(cfg)
    assert(cfg.K >= 2 && cfg.K == floor(cfg.K), ...
        "K must be an integer of at least two.");
    assert(numel(cfg.x) == cfg.K && all(isfinite(cfg.x)), ...
        "The spectral axis must contain K finite points.");
    assert(isfinite(cfg.gaussianWidth) && cfg.gaussianWidth > 0, ...
        "Gaussian width must be finite and positive.");
    assert(all(isfinite(cfg.separationValues)) && ...
        all(cfg.separationValues > 0), ...
        "Separations must be finite and positive.");
    assert(numel(cfg.trueAmplitude) == 2 && ...
        all(isfinite(cfg.trueAmplitude)) && all(cfg.trueAmplitude >= 0), ...
        "Two finite nonnegative true amplitudes are required.");
    assert(numel(cfg.echoTimesMs) >= 2 && ...
        all(isfinite(cfg.echoTimesMs)) && all(diff(cfg.echoTimesMs) > 0), ...
        "Echo times must be finite and strictly increasing.");
    assert(isfinite(cfg.sigmaNoise) && cfg.sigmaNoise >= 0, ...
        "Noise SD must be finite and nonnegative.");
    assert(cfg.nMonteCarlo >= 1 && ...
        cfg.nMonteCarlo == floor(cfg.nMonteCarlo), ...
        "Monte Carlo count must be a positive integer.");
    assert(size(cfg.trueT2Ms, 2) == 2 && ...
        all(isfinite(cfg.trueT2Ms), "all") && all(cfg.trueT2Ms > 0, "all"), ...
        "Each scenario must contain two finite positive T2 values.");
    assert(cfg.t2LowerBoundMs > 0 && ...
        cfg.t2UpperBoundMs > cfg.t2LowerBoundMs, ...
        "T2 bounds must define a strictly positive nonempty interval.");
end

function covariance = sampleCovariance(trajectory)
    centered = trajectory - mean(trajectory, 2);
    covariance = (centered*centered.')/(size(trajectory, 2) - 1);
    covariance = (covariance + covariance.')/2;
end

function ratio = theoreticalT2Ratio(t2Ms, echoTimesMs)
    decay = exp(-echoTimesMs/t2Ms);
    ratio = var(decay, 0)/(mean(decay)^2);
end

function fit = invertT2Ratio(observedRatio, echoTimesMs, lowerBound, ...
        upperBound, rootTolerance, boundTolerance)
    fit = struct("t2Ms", NaN, "valid", false, ...
        "boundHit", false, "statusCode", uint8(0));
    if ~isfinite(observedRatio)
        fit.statusCode = uint8(3);
        return;
    end
    if observedRatio < 0
        fit.statusCode = uint8(2);
        return;
    end

    ratioAtLower = theoreticalT2Ratio(lowerBound, echoTimesMs);
    ratioAtUpper = theoreticalT2Ratio(upperBound, echoTimesMs);
    theoreticalMinimum = min(ratioAtLower, ratioAtUpper);
    theoreticalMaximum = max(ratioAtLower, ratioAtUpper);
    if observedRatio < theoreticalMinimum || ...
            observedRatio > theoreticalMaximum
        fit.statusCode = uint8(4);
        return;
    end

    if observedRatio == ratioAtLower
        t2Ms = lowerBound;
        exitFlag = 1;
    elseif observedRatio == ratioAtUpper
        t2Ms = upperBound;
        exitFlag = 1;
    else
        rootFunction = @(candidateT2) theoreticalT2Ratio( ...
            candidateT2, echoTimesMs) - observedRatio;
        options = optimset("Display", "off", "TolX", rootTolerance);
        try
            [t2Ms, ~, exitFlag] = fzero(rootFunction, ...
                [lowerBound, upperBound], options);
        catch
            fit.statusCode = uint8(5);
            return;
        end
    end

    if exitFlag <= 0
        fit.statusCode = uint8(5);
        return;
    end
    if ~isfinite(t2Ms) || t2Ms <= 0
        fit.statusCode = uint8(3);
        return;
    end

    tolerance = boundTolerance*max(1, abs([lowerBound, upperBound]));
    fit.t2Ms = t2Ms;
    fit.valid = true;
    fit.boundHit = abs(t2Ms - lowerBound) <= tolerance(1) || ...
        abs(t2Ms - upperBound) <= tolerance(2);
    fit.statusCode = uint8(1);
end

function validation = runPreSimulationValidation(cfg, basisCache, ...
        trueTrajectories)
    nSeparations = numel(basisCache);
    nScenarios = size(trueTrajectories, 3);
    validation = struct();
    validation.projectionMaximumAbsoluteError = nan( ...
        nSeparations, nScenarios);
    validation.projectionPassed = false(nSeparations, nScenarios);
    validation.noiseCovarianceFrobeniusError = nan(nSeparations, 1);
    validation.noiseCovariancePassed = false(nSeparations, 1);
    validation.ratioRecoveredT2Ms = nan(nScenarios, 2);
    validation.ratioRecoveryAbsoluteErrorMs = nan(nScenarios, 2);
    validation.ratioRecoveryPassed = false(nScenarios, 2);

    for dIndex = 1:nSeparations
        cache = basisCache(dIndex);
        projectedNoiseCovariance = cfg.sigmaNoise^2* ...
            (cache.projection*cache.projection.');
        covarianceError = norm(projectedNoiseCovariance - ...
            cache.noiseCovariance, "fro");
        validation.noiseCovarianceFrobeniusError(dIndex) = ...
            covarianceError;
        validation.noiseCovariancePassed(dIndex) = covarianceError <= ...
            cfg.validation.noiseCovarianceTolerance;

        for scenarioIndex = 1:nScenarios
            amplitude = trueTrajectories(:, :, scenarioIndex);
            recovered = cache.projection*(cache.B*amplitude);
            projectionError = max(abs(recovered - amplitude), [], "all");
            validation.projectionMaximumAbsoluteError(dIndex, ...
                scenarioIndex) = projectionError;
            validation.projectionPassed(dIndex, scenarioIndex) = ...
                projectionError <= cfg.validation.projectionTolerance;
        end
    end

    for scenarioIndex = 1:nScenarios
        for metaboliteIndex = 1:2
            trajectory = trueTrajectories(metaboliteIndex, :, ...
                scenarioIndex).';
            observedRatio = var(trajectory, 0)/(mean(trajectory)^2);
            ratioFit = invertT2Ratio(observedRatio, cfg.echoTimesMs, ...
                cfg.t2LowerBoundMs, cfg.t2UpperBoundMs, ...
                cfg.rootToleranceMs, cfg.boundHitRelativeTolerance);
            validation.ratioRecoveredT2Ms(scenarioIndex, ...
                metaboliteIndex) = ratioFit.t2Ms;
            absoluteError = abs(ratioFit.t2Ms - ...
                cfg.trueT2Ms(scenarioIndex, metaboliteIndex));
            validation.ratioRecoveryAbsoluteErrorMs(scenarioIndex, ...
                metaboliteIndex) = absoluteError;
            validation.ratioRecoveryPassed(scenarioIndex, ...
                metaboliteIndex) = ratioFit.valid && absoluteError <= ...
                cfg.validation.ratioT2ToleranceMs;
        end
    end

    t2Grid = linspace(cfg.t2LowerBoundMs, cfg.t2UpperBoundMs, 2000);
    ratioGrid = arrayfun(@(t2) theoreticalT2Ratio( ...
        t2, cfg.echoTimesMs), t2Grid);
    validation.ratioMonotonicOnBounds = all(diff(ratioGrid) < 0);
end

function fit = solveBoundedLeastSquares(residualFunction, initialValue, ...
        lowerBound, upperBound, optimizerCfg, boundTolerance)
    initialValue = initialValue(:).';
    lowerBound = lowerBound(:).';
    upperBound = upperBound(:).';
    fit = failedNLSFit(numel(initialValue));

    try
        if exist("lsqnonlin", "file") == 2
            options = optimoptions("lsqnonlin", ...
                "Display", optimizerCfg.display, ...
                "FunctionTolerance", optimizerCfg.functionTolerance, ...
                "StepTolerance", optimizerCfg.stepTolerance, ...
                "OptimalityTolerance", optimizerCfg.optimalityTolerance, ...
                "MaxIterations", optimizerCfg.maxIterations, ...
                "MaxFunctionEvaluations", ...
                optimizerCfg.maxFunctionEvaluations);
            [candidate, ~, residual, exitFlag, output] = lsqnonlin( ...
                residualFunction, initialValue, lowerBound, upperBound, ...
                options);
            optimizerName = "lsqnonlin";
        else
            initialZ = boundedToUnconstrained( ...
                initialValue, lowerBound, upperBound);
            objective = @(z) sum(residualFunction( ...
                unconstrainedToBounded(z, lowerBound, upperBound)).^2);
            options = optimset("Display", char(optimizerCfg.display), ...
                "TolFun", optimizerCfg.functionTolerance, ...
                "TolX", optimizerCfg.stepTolerance, ...
                "MaxIter", optimizerCfg.maxIterations, ...
                "MaxFunEvals", optimizerCfg.maxFunctionEvaluations);
            [candidateZ, ~, exitFlag, output] = fminsearch( ...
                objective, initialZ, options);
            candidate = unconstrainedToBounded( ...
                candidateZ, lowerBound, upperBound);
            residual = residualFunction(candidate);
            optimizerName = "fminsearch-bounded-transform";
        end

        candidate = candidate(:).';
        fit.candidate = candidate;
        fit.exitFlag = exitFlag;
        fit.iterations = getOutputField(output, "iterations");
        fit.functionCount = getOutputField(output, "funcCount");
        fit.optimizerName = optimizerName;
        fit.boundHit = parameterAtBound( ...
            candidate, lowerBound, upperBound, boundTolerance);

        if ~all(isfinite(candidate)) || ~all(isfinite(residual))
            fit.statusCode = uint8(3);
        elseif exitFlag <= 0
            fit.statusCode = uint8(6);
        else
            fit.success = true;
            fit.statusCode = uint8(1);
        end
    catch
        fit.statusCode = uint8(5);
    end
end

function fit = failedNLSFit(nParameters)
    fit = struct("candidate", nan(1, nParameters), ...
        "success", false, "exitFlag", NaN, "iterations", NaN, ...
        "functionCount", NaN, "optimizerName", "", ...
        "boundHit", false(1, nParameters), "statusCode", uint8(5));
end

function value = getOutputField(output, fieldName)
    if isstruct(output) && isfield(output, fieldName) && ...
            isscalar(output.(fieldName))
        value = output.(fieldName);
    else
        value = NaN;
    end
end

function isAtBound = parameterAtBound(value, lowerBound, upperBound, ...
        relativeTolerance)
    lowerTolerance = relativeTolerance*max(1, abs(lowerBound));
    upperTolerance = relativeTolerance*max(1, abs(upperBound));
    isAtBound = abs(value - lowerBound) <= lowerTolerance | ...
        abs(value - upperBound) <= upperTolerance;
end

function z = boundedToUnconstrained(value, lowerBound, upperBound)
    fraction = (value - lowerBound)./(upperBound - lowerBound);
    fraction = min(max(fraction, sqrt(eps)), 1 - sqrt(eps));
    z = log(fraction./(1 - fraction));
end

function value = unconstrainedToBounded(z, lowerBound, upperBound)
    fraction = 1./(1 + exp(-z));
    value = lowerBound + (upperBound - lowerBound).*fraction;
end

function summary = summarizeT2Results(raw, cfg)
    nScenarios = size(cfg.trueT2Ms, 1);
    nSeparations = numel(cfg.separationValues);
    nMethods = numel(raw.t2MethodNames);
    nRows = nScenarios*nSeparations*nMethods;

    ScenarioIndex = nan(nRows, 1);
    Scenario = strings(nRows, 1);
    Separation = nan(nRows, 1);
    BasisCorrelation = nan(nRows, 1);
    GramConditionNumber = nan(nRows, 1);
    Method = strings(nRows, 1);
    TrueT2_1_ms = nan(nRows, 1);
    TrueT2_2_ms = nan(nRows, 1);
    BiasConditional_T2_1_ms = nan(nRows, 1);
    BiasConditional_T2_2_ms = nan(nRows, 1);
    RMSEConditional_T2_1_ms = nan(nRows, 1);
    RMSEConditional_T2_2_ms = nan(nRows, 1);
    MedianAbsErrorConditional_T2_1_ms = nan(nRows, 1);
    MedianAbsErrorConditional_T2_2_ms = nan(nRows, 1);
    ValidRatePercent_T2_1 = nan(nRows, 1);
    ValidRatePercent_T2_2 = nan(nRows, 1);

    row = 0;
    for scenarioIndex = 1:nScenarios
        for dIndex = 1:nSeparations
            for methodIndex = 1:nMethods
                row = row + 1;
                ScenarioIndex(row) = scenarioIndex;
                Scenario(row) = cfg.scenarioNames(scenarioIndex);
                Separation(row) = cfg.separationValues(dIndex);
                BasisCorrelation(row) = raw.basisCorrelation(dIndex);
                GramConditionNumber(row) = raw.gramConditionNumber(dIndex);
                Method(row) = raw.t2MethodNames(methodIndex);
                TrueT2_1_ms(row) = cfg.trueT2Ms(scenarioIndex, 1);
                TrueT2_2_ms(row) = cfg.trueT2Ms(scenarioIndex, 2);

                for metaboliteIndex = 1:2
                    estimates = reshape(raw.t2HatMs(scenarioIndex, ...
                        dIndex, :, methodIndex, metaboliteIndex), [], 1);
                    valid = reshape(raw.t2EstimateValid(scenarioIndex, ...
                        dIndex, :, methodIndex, metaboliteIndex), [], 1) & ...
                        isfinite(estimates);
                    [bias, rmse, medianAbsoluteError] = errorMetrics( ...
                        estimates(valid), ...
                        cfg.trueT2Ms(scenarioIndex, metaboliteIndex));
                    if metaboliteIndex == 1
                        BiasConditional_T2_1_ms(row) = bias;
                        RMSEConditional_T2_1_ms(row) = rmse;
                        MedianAbsErrorConditional_T2_1_ms(row) = ...
                            medianAbsoluteError;
                        ValidRatePercent_T2_1(row) = 100*mean(valid);
                    else
                        BiasConditional_T2_2_ms(row) = bias;
                        RMSEConditional_T2_2_ms(row) = rmse;
                        MedianAbsErrorConditional_T2_2_ms(row) = ...
                            medianAbsoluteError;
                        ValidRatePercent_T2_2(row) = 100*mean(valid);
                    end
                end
            end
        end
    end

    summary = table(ScenarioIndex, Scenario, Separation, ...
        BasisCorrelation, GramConditionNumber, Method, ...
        TrueT2_1_ms, TrueT2_2_ms, ...
        BiasConditional_T2_1_ms, BiasConditional_T2_2_ms, ...
        RMSEConditional_T2_1_ms, RMSEConditional_T2_2_ms, ...
        MedianAbsErrorConditional_T2_1_ms, ...
        MedianAbsErrorConditional_T2_2_ms, ...
        ValidRatePercent_T2_1, ValidRatePercent_T2_2);
end

function [summary, stats] = summarizeCovarianceResults(raw, cfg)
    nScenarios = size(cfg.trueT2Ms, 1);
    nSeparations = numel(cfg.separationValues);
    nMethods = numel(raw.covarianceMethodNames);
    nRows = nScenarios*nSeparations*nMethods;
    stats.meanEstimate = nan(2, 2, nScenarios, nSeparations, nMethods);
    stats.elementwiseBias = nan(2, 2, nScenarios, nSeparations, nMethods);
    stats.commonValidCount = nan(nScenarios, nSeparations);

    ScenarioIndex = nan(nRows, 1);
    Scenario = strings(nRows, 1);
    Separation = nan(nRows, 1);
    Method = strings(nRows, 1);
    CommonValidRepetitions = nan(nRows, 1);
    MeanS11 = nan(nRows, 1);
    MeanS12 = nan(nRows, 1);
    MeanS21 = nan(nRows, 1);
    MeanS22 = nan(nRows, 1);
    BiasS11 = nan(nRows, 1);
    BiasS12 = nan(nRows, 1);
    BiasS21 = nan(nRows, 1);
    BiasS22 = nan(nRows, 1);
    RMSE_Frobenius = nan(nRows, 1);

    row = 0;
    for scenarioIndex = 1:nScenarios
        reference = raw.trueCovariance(:, :, scenarioIndex);
        for dIndex = 1:nSeparations
            commonValid = true(cfg.nMonteCarlo, 1);
            for methodIndex = 1:nMethods
                estimates = raw.covarianceEstimate(:, :, scenarioIndex, ...
                    dIndex, :, methodIndex);
                methodValid = reshape(all(isfinite(estimates), [1, 2]), ...
                    [], 1);
                commonValid = commonValid & methodValid;
            end
            stats.commonValidCount(scenarioIndex, dIndex) = sum(commonValid);

            for methodIndex = 1:nMethods
                row = row + 1;
                ScenarioIndex(row) = scenarioIndex;
                Scenario(row) = cfg.scenarioNames(scenarioIndex);
                Separation(row) = cfg.separationValues(dIndex);
                Method(row) = raw.covarianceMethodNames(methodIndex);
                CommonValidRepetitions(row) = sum(commonValid);

                estimates = reshape(raw.covarianceEstimate(:, :, ...
                    scenarioIndex, dIndex, commonValid, methodIndex), ...
                    2, 2, []);
                if ~isempty(estimates)
                    meanEstimate = mean(estimates, 3);
                    elementwiseBias = meanEstimate - reference;
                    errors = estimates - reference;
                    squaredFrobenius = reshape(sum(errors.^2, [1, 2]), ...
                        [], 1);
                    rmseFrobenius = sqrt(mean(squaredFrobenius));
                    stats.meanEstimate(:, :, scenarioIndex, dIndex, ...
                        methodIndex) = meanEstimate;
                    stats.elementwiseBias(:, :, scenarioIndex, dIndex, ...
                        methodIndex) = elementwiseBias;
                    MeanS11(row) = meanEstimate(1, 1);
                    MeanS12(row) = meanEstimate(1, 2);
                    MeanS21(row) = meanEstimate(2, 1);
                    MeanS22(row) = meanEstimate(2, 2);
                    BiasS11(row) = elementwiseBias(1, 1);
                    BiasS12(row) = elementwiseBias(1, 2);
                    BiasS21(row) = elementwiseBias(2, 1);
                    BiasS22(row) = elementwiseBias(2, 2);
                    RMSE_Frobenius(row) = rmseFrobenius;
                end
            end
        end
    end

    summary = table(ScenarioIndex, Scenario, Separation, Method, ...
        CommonValidRepetitions, MeanS11, MeanS12, MeanS21, MeanS22, ...
        BiasS11, BiasS12, BiasS21, BiasS22, RMSE_Frobenius);
end

function summary = summarizeValidity(raw, cfg)
    nScenarios = size(cfg.trueT2Ms, 1);
    nSeparations = numel(cfg.separationValues);
    nMethods = numel(raw.t2MethodNames);
    nRows = nScenarios*nSeparations*nMethods*2;

    ScenarioIndex = nan(nRows, 1);
    Scenario = strings(nRows, 1);
    Separation = nan(nRows, 1);
    Method = strings(nRows, 1);
    Metabolite = nan(nRows, 1);
    Total = repmat(cfg.nMonteCarlo, nRows, 1);
    Valid = nan(nRows, 1);
    Invalid = nan(nRows, 1);
    ValidRatePercent = nan(nRows, 1);
    NegativeRatio = zeros(nRows, 1);
    Nonfinite = zeros(nRows, 1);
    OutsideTheoreticalRange = zeros(nRows, 1);
    RootOrSolverFailure = zeros(nRows, 1);
    OptimizerNonconvergence = zeros(nRows, 1);
    T2BoundHits = zeros(nRows, 1);
    AmplitudeBoundHits = nan(nRows, 1);

    row = 0;
    for scenarioIndex = 1:nScenarios
        for dIndex = 1:nSeparations
            for methodIndex = 1:nMethods
                for metaboliteIndex = 1:2
                    row = row + 1;
                    codes = reshape(raw.t2StatusCode(scenarioIndex, ...
                        dIndex, :, methodIndex, metaboliteIndex), [], 1);
                    valid = reshape(raw.t2EstimateValid(scenarioIndex, ...
                        dIndex, :, methodIndex, metaboliteIndex), [], 1);
                    boundHits = reshape(raw.t2BoundHit(scenarioIndex, ...
                        dIndex, :, methodIndex, metaboliteIndex), [], 1);
                    ScenarioIndex(row) = scenarioIndex;
                    Scenario(row) = cfg.scenarioNames(scenarioIndex);
                    Separation(row) = cfg.separationValues(dIndex);
                    Method(row) = raw.t2MethodNames(methodIndex);
                    Metabolite(row) = metaboliteIndex;
                    Valid(row) = sum(valid);
                    Invalid(row) = cfg.nMonteCarlo - Valid(row);
                    ValidRatePercent(row) = 100*mean(valid);
                    NegativeRatio(row) = sum(codes == 2);
                    Nonfinite(row) = sum(codes == 3);
                    OutsideTheoreticalRange(row) = sum(codes == 4);
                    RootOrSolverFailure(row) = sum(codes == 5);
                    OptimizerNonconvergence(row) = sum(codes == 6);
                    T2BoundHits(row) = sum(boundHits(valid) == 1);
                    if methodIndex == 3
                        amplitudeBound = reshape( ...
                            raw.nlsAmplitudeBoundHit(scenarioIndex, ...
                            dIndex, :, metaboliteIndex), [], 1);
                        AmplitudeBoundHits(row) = sum(amplitudeBound(valid));
                    end
                end
            end
        end
    end

    summary = table(ScenarioIndex, Scenario, Separation, Method, ...
        Metabolite, Total, Valid, Invalid, ValidRatePercent, ...
        NegativeRatio, Nonfinite, OutsideTheoreticalRange, ...
        RootOrSolverFailure, OptimizerNonconvergence, T2BoundHits, ...
        AmplitudeBoundHits);
end

function summary = summarizeCorrectedCovarianceDiagnostics(raw, cfg)
    nScenarios = size(cfg.trueT2Ms, 1);
    nSeparations = numel(cfg.separationValues);
    nRows = nScenarios*nSeparations;
    ScenarioIndex = nan(nRows, 1);
    Scenario = strings(nRows, 1);
    Separation = nan(nRows, 1);
    NegativeDiagonalCount = nan(nRows, 1);
    IndefiniteCount = nan(nRows, 1);
    MinimumObservedEigenvalue = nan(nRows, 1);

    row = 0;
    for scenarioIndex = 1:nScenarios
        for dIndex = 1:nSeparations
            row = row + 1;
            negativeDiagonal = reshape( ...
                raw.correctedHasNegativeDiagonal(scenarioIndex, ...
                dIndex, :), [], 1);
            minimumEigenvalue = reshape(raw.correctedMinEigenvalue( ...
                scenarioIndex, dIndex, :), [], 1);
            ScenarioIndex(row) = scenarioIndex;
            Scenario(row) = cfg.scenarioNames(scenarioIndex);
            Separation(row) = cfg.separationValues(dIndex);
            NegativeDiagonalCount(row) = sum(negativeDiagonal);
            IndefiniteCount(row) = sum(minimumEigenvalue < 0);
            MinimumObservedEigenvalue(row) = min(minimumEigenvalue);
        end
    end
    summary = table(ScenarioIndex, Scenario, Separation, ...
        NegativeDiagonalCount, IndefiniteCount, ...
        MinimumObservedEigenvalue);
end

function summary = validateExpectedCovarianceBias(raw, cfg, basisCache)
    nScenarios = size(cfg.trueT2Ms, 1);
    nSeparations = numel(cfg.separationValues);
    nRows = nScenarios*nSeparations*2;
    ScenarioIndex = nan(nRows, 1);
    Scenario = strings(nRows, 1);
    Separation = nan(nRows, 1);
    Method = strings(nRows, 1);
    MaximumAbsoluteDeviation = nan(nRows, 1);
    MaximumStandardizedDeviation = nan(nRows, 1);
    StandardErrorMultiplier = repmat( ...
        cfg.validation.standardErrorMultiplier, nRows, 1);
    Evaluated = false(nRows, 1);
    Passed = false(nRows, 1);

    row = 0;
    for scenarioIndex = 1:nScenarios
        for dIndex = 1:nSeparations
            for methodIndex = 1:2
                row = row + 1;
                estimates = reshape(raw.covarianceEstimate(:, :, ...
                    scenarioIndex, dIndex, :, methodIndex), 2, 2, []);
                if methodIndex == 1
                    expected = raw.trueCovariance(:, :, scenarioIndex) + ...
                        basisCache(dIndex).noiseCovariance;
                else
                    expected = raw.trueCovariance(:, :, scenarioIndex);
                end
                meanEstimate = mean(estimates, 3);
                deviation = meanEstimate - expected;
                standardError = std(estimates, 0, 3)/sqrt(cfg.nMonteCarlo);
                standardized = abs(deviation)./standardError;
                bothNumericallyZero = standardError == 0 & ...
                    abs(deviation) <= cfg.validation.noiseCovarianceTolerance;
                standardized(bothNumericallyZero) = 0;

                ScenarioIndex(row) = scenarioIndex;
                Scenario(row) = cfg.scenarioNames(scenarioIndex);
                Separation(row) = cfg.separationValues(dIndex);
                Method(row) = raw.covarianceMethodNames(methodIndex);
                MaximumAbsoluteDeviation(row) = max(abs(deviation), [], "all");
                MaximumStandardizedDeviation(row) = ...
                    max(standardized, [], "all");
                Evaluated(row) = cfg.nMonteCarlo >= ...
                    cfg.validation.minimumMonteCarloForSamplingCheck;
                Passed(row) = Evaluated(row) && ...
                    MaximumStandardizedDeviation(row) <= ...
                    cfg.validation.standardErrorMultiplier;
            end
        end
    end

    summary = table(ScenarioIndex, Scenario, Separation, Method, ...
        MaximumAbsoluteDeviation, MaximumStandardizedDeviation, ...
        StandardErrorMultiplier, Evaluated, Passed);
end

function [bias, rmse, medianAbsoluteError] = errorMetrics(estimates, truth)
    if isempty(estimates)
        bias = NaN;
        rmse = NaN;
        medianAbsoluteError = NaN;
        return;
    end
    errors = estimates - truth;
    bias = mean(errors);
    rmse = sqrt(mean(errors.^2));
    medianAbsoluteError = median(abs(errors));
end

function saveFigureFiles(fig, cfg, baseName)
    if ~cfg.saveOutputs
        return;
    end
    exportgraphics(fig, fullfile(cfg.outputDirectory, baseName + ".png"), ...
        "Resolution", 200);
    savefig(fig, fullfile(cfg.outputDirectory, baseName + ".fig"));
end

function addMatrixValueLabels(ax, matrix, colorLimits)
    midpoint = mean(colorLimits);
    for row = 1:2
        for column = 1:2
            if matrix(row, column) < midpoint
                textColor = "w";
            else
                textColor = "k";
            end
            text(ax, column, row, sprintf("%.4g", matrix(row, column)), ...
                "HorizontalAlignment", "center", ...
                "VerticalAlignment", "middle", ...
                "Color", textColor, "FontWeight", "bold", ...
                "FontSize", 10);
        end
    end
end

function stylePlotAxes(ax)
    set(ax, "Color", "w", "XColor", "k", "YColor", "k", ...
        "GridColor", [0.65, 0.65, 0.65], "GridAlpha", 0.45, ...
        "FontSize", 11);
    ax.Title.Color = "k";
    ax.XLabel.Color = "k";
    ax.YLabel.Color = "k";
end
