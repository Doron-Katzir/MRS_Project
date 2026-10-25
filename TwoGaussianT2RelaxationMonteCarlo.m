%% TwoGaussianT2RelaxationMonteCarlo
% Standalone Monte Carlo comparison of three T2 estimators for two
% overlapping Gaussian spectral basis functions.
%
% The simulated K-by-N data matrix is
%
%   S(:,i) = A1*exp(-TE(i)/T2_1)*b1 ...
%          + A2*exp(-TE(i)/T2_2)*b2 + epsilon(:,i),
%
% with independent epsilon(:,i) ~ N(0, sigmaNoise^2*I_K). There are no
% random effects or concentration fluctuations in this simulation.

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

% Fixed guesses are deliberately different from every true parameter set.
% The same corresponding guesses and bounds are used in the naive and
% joint nonlinear fits.
cfg.initialAmplitude = [0.8, 0.8];
cfg.initialT2Ms = [180, 180];
cfg.amplitudeLowerBound = [0, 0];
cfg.amplitudeUpperBound = [5, 5];
cfg.t2LowerBoundMs = [1, 1];
cfg.t2UpperBoundMs = [2000, 2000];

% Both nonlinear methods receive these same stopping tolerances.
cfg.optimizer.functionTolerance = 1e-10;
cfg.optimizer.stepTolerance = 1e-10;
cfg.optimizer.optimalityTolerance = 1e-8;
cfg.optimizer.maxIterations = 500;
cfg.optimizer.maxFunctionEvaluations = 4000;
cfg.optimizer.display = "off";
cfg.boundHitRelativeTolerance = 1e-6;

cfg.makeFigures = true;
cfg.saveOutputs = true;
cfg.figureVisible = "on";
cfg.outputDirectory = fullfile(fileparts(mfilename("fullpath")), ...
    "T2RelaxationMonteCarloResults");

% Set MRS_T2_SMOKE_TEST=1 only for a fast automated code check. It changes
% the run to two repetitions and two separations, and suppresses outputs.
if strcmp(getenv("MRS_T2_SMOKE_TEST"), "1")
    cfg.nMonteCarlo = 2;
    cfg.separationValues = cfg.separationValues(1:2);
    cfg.makeFigures = false;
    cfg.saveOutputs = false;
end

validateConfiguration(cfg);
rng(cfg.rngSeed, "twister");

methodNames = ["Naive OLS", "Corrected OLS", "Joint NLS"];
nMethods = numel(methodNames);
nScenarios = size(cfg.trueT2Ms, 1);
nSeparations = numel(cfg.separationValues);
nEchoes = numel(cfg.echoTimesMs);

%% Precompute bases, OLS projections, and known OLS noise variances

basisCache = repmat(struct( ...
    "B", [], "projection", [], "amplitudeNoiseVariance", [], ...
    "basisCorrelation", NaN, "gramConditionNumber", NaN), ...
    nSeparations, 1);

for dIndex = 1:nSeparations
    d = cfg.separationValues(dIndex);
    centers = [5 - d/2, 5 + d/2];
    b1 = exp(-0.5*((cfg.x - centers(1))/cfg.gaussianWidth).^2);
    b2 = exp(-0.5*((cfg.x - centers(2))/cfg.gaussianWidth).^2);
    b1 = b1/max(b1);  % exact unit peak height on the sampled grid
    b2 = b2/max(b2);
    B = [b1, b2];
    gram = B.'*B;

    % Linear solves are used throughout; no explicit matrix inverse is
    % formed. B\eye(K) is the precomputed two-channel OLS projection.
    basisCache(dIndex).B = B;
    basisCache(dIndex).projection = B\eye(cfg.K);
    basisCache(dIndex).amplitudeNoiseVariance = ...
        cfg.sigmaNoise^2*diag(gram\eye(2));
    basisCache(dIndex).basisCorrelation = ...
        (b1.'*b2)/sqrt((b1.'*b1)*(b2.'*b2));
    basisCache(dIndex).gramConditionNumber = cond(gram);
end

%% Allocate raw Monte Carlo results

resultSize = [nScenarios, nSeparations, cfg.nMonteCarlo, ...
    nMethods, 2];

raw = struct();
raw.methodNames = methodNames;
raw.scenarioNames = cfg.scenarioNames;
raw.separationValues = cfg.separationValues;
raw.basisCorrelation = reshape( ...
    [basisCache.basisCorrelation], [], 1);
raw.gramConditionNumber = reshape( ...
    [basisCache.gramConditionNumber], [], 1);
raw.t2HatMs = nan(resultSize);
raw.amplitudeHat = nan(resultSize);
raw.candidateT2Ms = nan(resultSize);
raw.candidateAmplitude = nan(resultSize);
raw.estimateValid = false(resultSize);

% Values are NaN for the corrected moment estimator because it does not
% invoke an optimizer or have optimization bounds.
raw.optimizerConverged = nan(resultSize);
raw.exitFlag = nan(resultSize);
raw.iterations = nan(resultSize);
raw.functionCount = nan(resultSize);
raw.t2BoundHit = nan(resultSize);
raw.amplitudeBoundHit = nan(resultSize);

% Corrected-status codes: 0 not evaluated, 1 valid, 2 denominator <= 0,
% 3 nonfinite q, and 4 q outside (0,1).
raw.correctedStatusCode = zeros( ...
    nScenarios, nSeparations, cfg.nMonteCarlo, 2, "uint8");
raw.correctedQHat = nan( ...
    nScenarios, nSeparations, cfg.nMonteCarlo, 2);
raw.correctedDenominator = nan( ...
    nScenarios, nSeparations, cfg.nMonteCarlo, 2);

% Nonlinear failure codes: 0 success, 1 nonpositive exit flag,
% 2 nonfinite solution/residual, and 3 solver exception.
raw.nonlinearFailureCode = zeros(resultSize, "uint8");
raw.optimizerName = strings( ...
    nScenarios, nSeparations, cfg.nMonteCarlo, nMethods);

%% Monte Carlo experiment

deltaTE = cfg.echoTimesMs(2) - cfg.echoTimesMs(1);
initialTheta = [cfg.initialAmplitude, cfg.initialT2Ms];
jointLower = [cfg.amplitudeLowerBound, cfg.t2LowerBoundMs];
jointUpper = [cfg.amplitudeUpperBound, cfg.t2UpperBoundMs];

for scenarioIndex = 1:nScenarios
    for dIndex = 1:nSeparations
        cache = basisCache(dIndex);
        trueTrajectories = cfg.trueAmplitude(:).*exp( ...
            -cfg.echoTimesMs.'./cfg.trueT2Ms(scenarioIndex, :).');
        cleanSpectra = cache.B*trueTrajectories;

        fprintf("Scenario %d/%d, d = %.3f (%d/%d), cond(B''B) = %.3g\n", ...
            scenarioIndex, nScenarios, cfg.separationValues(dIndex), ...
            dIndex, nSeparations, cache.gramConditionNumber);

        for mcIndex = 1:cfg.nMonteCarlo
            % This single noisy dataset is reused by all three methods.
            spectra = cleanSpectra + ...
                cfg.sigmaNoise*randn(cfg.K, nEchoes);
            amplitudeOLS = cache.projection*spectra;

            % Method A: independent nonlinear fits to the two OLS
            % amplitude trajectories, without logarithmic transformation.
            for metaboliteIndex = 1:2
                residualFunction = @(p) p(1)*exp( ...
                    -cfg.echoTimesMs/p(2)) - ...
                    amplitudeOLS(metaboliteIndex, :).';
                fit = solveBoundedLeastSquares(residualFunction, ...
                    [cfg.initialAmplitude(metaboliteIndex), ...
                    cfg.initialT2Ms(metaboliteIndex)], ...
                    [cfg.amplitudeLowerBound(metaboliteIndex), ...
                    cfg.t2LowerBoundMs(metaboliteIndex)], ...
                    [cfg.amplitudeUpperBound(metaboliteIndex), ...
                    cfg.t2UpperBoundMs(metaboliteIndex)], cfg.optimizer, ...
                    cfg.boundHitRelativeTolerance);

                raw = storeNonlinearFit(raw, fit, scenarioIndex, ...
                    dIndex, mcIndex, 1, metaboliteIndex, 1, 2);
            end

            % Method B: corrected adjacent-echo moment estimator, using
            % the known simulation noise variance propagated through OLS.
            for metaboliteIndex = 1:2
                trajectory = amplitudeOLS(metaboliteIndex, :).';
                corrected = correctedMomentT2(trajectory, deltaTE, ...
                    cache.amplitudeNoiseVariance(metaboliteIndex));
                raw.correctedStatusCode(scenarioIndex, dIndex, ...
                    mcIndex, metaboliteIndex) = corrected.statusCode;
                raw.correctedQHat(scenarioIndex, dIndex, ...
                    mcIndex, metaboliteIndex) = corrected.qHat;
                raw.correctedDenominator(scenarioIndex, dIndex, ...
                    mcIndex, metaboliteIndex) = corrected.denominator;
                if corrected.valid
                    raw.t2HatMs(scenarioIndex, dIndex, mcIndex, 2, ...
                        metaboliteIndex) = corrected.t2Ms;
                    raw.estimateValid(scenarioIndex, dIndex, mcIndex, 2, ...
                        metaboliteIndex) = true;
                end
            end

            % Method C: one simultaneous fit to every spectral point and
            % echo. Its residual vector contains K*N entries.
            jointResidual = @(theta) jointSpectralResidual(theta, ...
                cfg.echoTimesMs, cache.B, spectra);
            jointFit = solveBoundedLeastSquares(jointResidual, ...
                initialTheta, jointLower, jointUpper, cfg.optimizer, ...
                cfg.boundHitRelativeTolerance);
            for metaboliteIndex = 1:2
                raw = storeNonlinearFit(raw, jointFit, scenarioIndex, ...
                    dIndex, mcIndex, 3, metaboliteIndex, ...
                    metaboliteIndex, metaboliteIndex + 2);
            end
        end
    end
end

%% Summaries and failure/bound reports

summaryTable = summarizeResults(raw, cfg);
statusTable = summarizeCorrectedStatus(raw, cfg);

disp(summaryTable);
disp("Corrected-estimator validity status counts:");
disp(statusTable);

if cfg.saveOutputs
    if ~isfolder(cfg.outputDirectory)
        mkdir(cfg.outputDirectory);
    end
    save(fullfile(cfg.outputDirectory, "TwoGaussianT2_RawResults.mat"), ...
        "raw", "summaryTable", "statusTable", "cfg", "-v7.3");
    writetable(summaryTable, fullfile(cfg.outputDirectory, ...
        "TwoGaussianT2_Summary.csv"));
    writetable(statusTable, fullfile(cfg.outputDirectory, ...
        "TwoGaussianT2_CorrectedStatus.csv"));
end

%% Accuracy figures

if cfg.makeFigures
    colors = lines(nMethods);
    markers = ["o", "s", "^"];

    for scenarioIndex = 1:nScenarios
        fig = figure("Color", "w", "Visible", cfg.figureVisible, ...
            "Name", cfg.scenarioNames(scenarioIndex));
        layout = tiledlayout(fig, 2, 2, "TileSpacing", "compact", ...
            "Padding", "compact");

        for metaboliteIndex = 1:2
            ax = nexttile(layout, metaboliteIndex);
            hold(ax, "on");
            for methodIndex = 1:nMethods
                rows = summaryTable.ScenarioIndex == scenarioIndex & ...
                    summaryTable.Method == methodNames(methodIndex);
                if metaboliteIndex == 1
                    y = summaryTable.RMSE_T2_1_ms(rows);
                else
                    y = summaryTable.RMSE_T2_2_ms(rows);
                end
                plot(ax, summaryTable.Separation(rows), y, ...
                    "-" + markers(methodIndex), ...
                    "Color", colors(methodIndex, :), ...
                    "LineWidth", 1.7, "MarkerSize", 6, ...
                    "DisplayName", methodNames(methodIndex));
            end
            xlabel(ax, "Gaussian center separation d");
            ylabel(ax, sprintf("T2_%d RMSE (ms)", metaboliteIndex));
            title(ax, sprintf("T2_%d = %.0f ms", metaboliteIndex, ...
                cfg.trueT2Ms(scenarioIndex, metaboliteIndex)));
            grid(ax, "on");
            if metaboliteIndex == 1
                legend(ax, "Location", "best");
            end

            validityAx = nexttile(layout, 2 + metaboliteIndex);
            hold(validityAx, "on");
            for methodIndex = 1:nMethods
                rows = summaryTable.ScenarioIndex == scenarioIndex & ...
                    summaryTable.Method == methodNames(methodIndex);
                if metaboliteIndex == 1
                    validRate = summaryTable.ValidRate_T2_1(rows);
                else
                    validRate = summaryTable.ValidRate_T2_2(rows);
                end
                plot(validityAx, summaryTable.Separation(rows), validRate, ...
                    "-" + markers(methodIndex), ...
                    "Color", colors(methodIndex, :), ...
                    "LineWidth", 1.7, "MarkerSize", 6, ...
                    "DisplayName", methodNames(methodIndex));
            end
            xlabel(validityAx, "Gaussian center separation d");
            ylabel(validityAx, sprintf("Valid T2_%d fraction", ...
                metaboliteIndex));
            ylim(validityAx, [0, 1.05]);
            grid(validityAx, "on");
        end
        title(layout, cfg.scenarioNames(scenarioIndex));

        if cfg.saveOutputs
            baseName = sprintf("TwoGaussianT2_RMSE_Scenario%d", ...
                scenarioIndex);
            exportgraphics(fig, fullfile(cfg.outputDirectory, ...
                baseName + ".png"), "Resolution", 200);
            savefig(fig, fullfile(cfg.outputDirectory, baseName + ".fig"));
        end
    end
end

%% Local functions

function validateConfiguration(cfg)
    assert(cfg.K >= 2 && cfg.K == floor(cfg.K), ...
        "K must be an integer of at least two.");
    assert(numel(cfg.x) == cfg.K && all(isfinite(cfg.x)), ...
        "The spectral axis must contain K finite points.");
    assert(cfg.gaussianWidth > 0 && isfinite(cfg.gaussianWidth), ...
        "Gaussian width must be finite and positive.");
    assert(all(isfinite(cfg.separationValues)) && ...
        all(cfg.separationValues > 0), ...
        "All Gaussian separations must be finite and positive.");
    assert(numel(cfg.echoTimesMs) >= 2 && ...
        all(diff(cfg.echoTimesMs) > 0), ...
        "Echo times must be finite and strictly increasing.");
    spacing = diff(cfg.echoTimesMs);
    spacingTolerance = 10*eps(max(abs(cfg.echoTimesMs)));
    assert(all(abs(spacing - spacing(1)) <= spacingTolerance), ...
        "The corrected estimator requires equally spaced echo times.");
    assert(cfg.sigmaNoise >= 0 && isfinite(cfg.sigmaNoise), ...
        "Noise standard deviation must be finite and nonnegative.");
    assert(cfg.nMonteCarlo >= 1 && ...
        cfg.nMonteCarlo == floor(cfg.nMonteCarlo), ...
        "nMonteCarlo must be a positive integer.");
    assert(size(cfg.trueT2Ms, 2) == 2 && ...
        all(isfinite(cfg.trueT2Ms), "all") && all(cfg.trueT2Ms > 0, "all"), ...
        "Each scenario must contain two finite positive T2 values.");
end

function residual = jointSpectralResidual(theta, echoTimesMs, B, spectra)
    amplitude = theta(1:2).';
    t2Ms = theta(3:4).';
    trajectories = amplitude.*exp(-echoTimesMs.'./t2Ms);
    modelSpectra = B*trajectories;
    residual = modelSpectra(:) - spectra(:);
end

function corrected = correctedMomentT2(amplitudeTrajectory, deltaTE, ...
        amplitudeNoiseVariance)
    nEchoes = numel(amplitudeTrajectory);
    numerator = sum(amplitudeTrajectory(1:end-1).* ...
        amplitudeTrajectory(2:end));
    denominator = sum(amplitudeTrajectory(1:end-1).^2) - ...
        (nEchoes - 1)*amplitudeNoiseVariance;

    corrected = struct("t2Ms", NaN, "qHat", NaN, ...
        "denominator", denominator, "valid", false, ...
        "statusCode", uint8(0));
    if ~isfinite(denominator) || denominator <= 0
        corrected.statusCode = uint8(2);
        return;
    end

    qHat = numerator/denominator;
    corrected.qHat = qHat;
    if ~isfinite(qHat)
        corrected.statusCode = uint8(3);
        return;
    end
    if qHat <= 0 || qHat >= 1
        corrected.statusCode = uint8(4);
        return;
    end

    t2Ms = -deltaTE/log(qHat);
    if ~isfinite(t2Ms) || t2Ms <= 0
        corrected.statusCode = uint8(3);
        return;
    end
    corrected.t2Ms = t2Ms;
    corrected.valid = true;
    corrected.statusCode = uint8(1);
end

function fit = solveBoundedLeastSquares(residualFunction, initialValue, ...
        lowerBound, upperBound, optimizerCfg, boundTolerance)
    initialValue = initialValue(:).';
    lowerBound = lowerBound(:).';
    upperBound = upperBound(:).';
    fit = failedFit(numel(initialValue), uint8(3));

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
            % Base-MATLAB fallback: a logistic transformation enforces the
            % same finite bounds, and fminsearch minimizes residual sum of
            % squares. Both naive and joint fits use this same wrapper.
            initialZ = boundedToUnconstrained(initialValue, ...
                lowerBound, upperBound);
            objective = @(z) sum(residualFunction( ...
                unconstrainedToBounded(z, lowerBound, upperBound)).^2);
            options = optimset("Display", char(optimizerCfg.display), ...
                "TolFun", optimizerCfg.functionTolerance, ...
                "TolX", optimizerCfg.stepTolerance, ...
                "MaxIter", optimizerCfg.maxIterations, ...
                "MaxFunEvals", optimizerCfg.maxFunctionEvaluations);
            [candidateZ, ~, exitFlag, output] = fminsearch( ...
                objective, initialZ, options);
            candidate = unconstrainedToBounded(candidateZ, ...
                lowerBound, upperBound);
            residual = residualFunction(candidate);
            optimizerName = "fminsearch-bounded-transform";
        end

        candidate = candidate(:).';
        residualIsFinite = all(isfinite(residual));
        candidateIsFinite = all(isfinite(candidate));
        fit.candidate = candidate;
        fit.exitFlag = exitFlag;
        fit.iterations = getOutputField(output, "iterations");
        fit.functionCount = getOutputField(output, "funcCount");
        fit.optimizerName = optimizerName;
        fit.boundHit = parameterAtBound(candidate, lowerBound, ...
            upperBound, boundTolerance);

        if ~candidateIsFinite || ~residualIsFinite
            fit.failureCode = uint8(2);
        elseif exitFlag <= 0
            fit.failureCode = uint8(1);
        else
            fit.success = true;
            fit.failureCode = uint8(0);
        end
    catch
        % The caller retains NaNs and failure code 3 for this repetition.
    end
end

function value = getOutputField(output, fieldName)
    if isstruct(output) && isfield(output, fieldName) && ...
            isscalar(output.(fieldName))
        value = output.(fieldName);
    else
        value = NaN;
    end
end

function fit = failedFit(nParameters, failureCode)
    fit = struct("candidate", nan(1, nParameters), ...
        "success", false, "exitFlag", NaN, "iterations", NaN, ...
        "functionCount", NaN, "optimizerName", "", ...
        "boundHit", false(1, nParameters), ...
        "failureCode", failureCode);
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

function raw = storeNonlinearFit(raw, fit, scenarioIndex, dIndex, ...
        mcIndex, methodIndex, metaboliteIndex, amplitudeParameterIndex, ...
        t2ParameterIndex)
    index = {scenarioIndex, dIndex, mcIndex, methodIndex, ...
        metaboliteIndex};
    raw.candidateAmplitude(index{:}) = ...
        fit.candidate(amplitudeParameterIndex);
    raw.candidateT2Ms(index{:}) = fit.candidate(t2ParameterIndex);
    raw.optimizerConverged(index{:}) = double(fit.success);
    raw.exitFlag(index{:}) = fit.exitFlag;
    raw.iterations(index{:}) = fit.iterations;
    raw.functionCount(index{:}) = fit.functionCount;
    raw.amplitudeBoundHit(index{:}) = double( ...
        fit.boundHit(amplitudeParameterIndex));
    raw.t2BoundHit(index{:}) = double(fit.boundHit(t2ParameterIndex));
    raw.nonlinearFailureCode(index{:}) = fit.failureCode;
    raw.optimizerName(scenarioIndex, dIndex, mcIndex, methodIndex) = ...
        fit.optimizerName;

    if fit.success
        raw.amplitudeHat(index{:}) = fit.candidate(amplitudeParameterIndex);
        raw.t2HatMs(index{:}) = fit.candidate(t2ParameterIndex);
        raw.estimateValid(index{:}) = true;
    end
end

function summary = summarizeResults(raw, cfg)
    nScenarios = size(cfg.trueT2Ms, 1);
    nSeparations = numel(cfg.separationValues);
    nMethods = numel(raw.methodNames);
    nRows = nScenarios*nSeparations*nMethods;

    ScenarioIndex = nan(nRows, 1);
    Scenario = strings(nRows, 1);
    Separation = nan(nRows, 1);
    BasisCorrelation = nan(nRows, 1);
    GramConditionNumber = nan(nRows, 1);
    Method = strings(nRows, 1);
    TrueT2_1_ms = nan(nRows, 1);
    TrueT2_2_ms = nan(nRows, 1);
    Bias_T2_1_ms = nan(nRows, 1);
    Bias_T2_2_ms = nan(nRows, 1);
    RMSE_T2_1_ms = nan(nRows, 1);
    RMSE_T2_2_ms = nan(nRows, 1);
    MedianAbsError_T2_1_ms = nan(nRows, 1);
    MedianAbsError_T2_2_ms = nan(nRows, 1);
    ValidRate_T2_1 = nan(nRows, 1);
    ValidRate_T2_2 = nan(nRows, 1);
    BothT2ValidRate = nan(nRows, 1);
    OptimizerFailureRate_T2_1 = nan(nRows, 1);
    OptimizerFailureRate_T2_2 = nan(nRows, 1);
    AmplitudeBoundHitRate_1 = nan(nRows, 1);
    AmplitudeBoundHitRate_2 = nan(nRows, 1);
    T2BoundHitRate_1 = nan(nRows, 1);
    T2BoundHitRate_2 = nan(nRows, 1);

    row = 0;
    for scenarioIndex = 1:nScenarios
        for dIndex = 1:nSeparations
            for methodIndex = 1:nMethods
                row = row + 1;
                ScenarioIndex(row) = scenarioIndex;
                Scenario(row) = cfg.scenarioNames(scenarioIndex);
                Separation(row) = cfg.separationValues(dIndex);
                BasisCorrelation(row) = raw.basisCorrelation(dIndex);
                GramConditionNumber(row) = ...
                    raw.gramConditionNumber(dIndex);
                Method(row) = raw.methodNames(methodIndex);
                TrueT2_1_ms(row) = cfg.trueT2Ms(scenarioIndex, 1);
                TrueT2_2_ms(row) = cfg.trueT2Ms(scenarioIndex, 2);

                estimate1 = reshape(raw.t2HatMs(scenarioIndex, ...
                    dIndex, :, methodIndex, 1), [], 1);
                estimate2 = reshape(raw.t2HatMs(scenarioIndex, ...
                    dIndex, :, methodIndex, 2), [], 1);
                valid1 = isfinite(estimate1);
                valid2 = isfinite(estimate2);
                ValidRate_T2_1(row) = mean(valid1);
                ValidRate_T2_2(row) = mean(valid2);
                BothT2ValidRate(row) = mean(valid1 & valid2);

                [Bias_T2_1_ms(row), RMSE_T2_1_ms(row), ...
                    MedianAbsError_T2_1_ms(row)] = errorMetrics( ...
                    estimate1(valid1), TrueT2_1_ms(row));
                [Bias_T2_2_ms(row), RMSE_T2_2_ms(row), ...
                    MedianAbsError_T2_2_ms(row)] = errorMetrics( ...
                    estimate2(valid2), TrueT2_2_ms(row));

                if methodIndex ~= 2
                    converged1 = reshape(raw.optimizerConverged( ...
                        scenarioIndex, dIndex, :, methodIndex, 1), [], 1);
                    converged2 = reshape(raw.optimizerConverged( ...
                        scenarioIndex, dIndex, :, methodIndex, 2), [], 1);
                    OptimizerFailureRate_T2_1(row) = ...
                        mean(converged1 == 0);
                    OptimizerFailureRate_T2_2(row) = ...
                        mean(converged2 == 0);
                    AmplitudeBoundHitRate_1(row) = meanForMask( ...
                        raw.amplitudeBoundHit(scenarioIndex, dIndex, :, ...
                        methodIndex, 1), converged1 == 1);
                    AmplitudeBoundHitRate_2(row) = meanForMask( ...
                        raw.amplitudeBoundHit(scenarioIndex, dIndex, :, ...
                        methodIndex, 2), converged2 == 1);
                    T2BoundHitRate_1(row) = meanForMask( ...
                        raw.t2BoundHit(scenarioIndex, dIndex, :, ...
                        methodIndex, 1), converged1 == 1);
                    T2BoundHitRate_2(row) = meanForMask( ...
                        raw.t2BoundHit(scenarioIndex, dIndex, :, ...
                        methodIndex, 2), converged2 == 1);
                end
            end
        end
    end

    summary = table(ScenarioIndex, Scenario, Separation, ...
        BasisCorrelation, GramConditionNumber, Method, ...
        TrueT2_1_ms, TrueT2_2_ms, Bias_T2_1_ms, Bias_T2_2_ms, ...
        RMSE_T2_1_ms, RMSE_T2_2_ms, ...
        MedianAbsError_T2_1_ms, MedianAbsError_T2_2_ms, ...
        ValidRate_T2_1, ValidRate_T2_2, BothT2ValidRate, ...
        OptimizerFailureRate_T2_1, OptimizerFailureRate_T2_2, ...
        AmplitudeBoundHitRate_1, AmplitudeBoundHitRate_2, ...
        T2BoundHitRate_1, T2BoundHitRate_2);
end

function status = summarizeCorrectedStatus(raw, cfg)
    nScenarios = size(cfg.trueT2Ms, 1);
    nSeparations = numel(cfg.separationValues);
    nRows = nScenarios*nSeparations*2;
    ScenarioIndex = nan(nRows, 1);
    Scenario = strings(nRows, 1);
    Separation = nan(nRows, 1);
    Metabolite = nan(nRows, 1);
    Valid = nan(nRows, 1);
    NonpositiveDenominator = nan(nRows, 1);
    NonfiniteQ = nan(nRows, 1);
    QOutsideOpenUnitInterval = nan(nRows, 1);

    row = 0;
    for scenarioIndex = 1:nScenarios
        for dIndex = 1:nSeparations
            for metaboliteIndex = 1:2
                row = row + 1;
                codes = reshape(raw.correctedStatusCode(scenarioIndex, ...
                    dIndex, :, metaboliteIndex), [], 1);
                ScenarioIndex(row) = scenarioIndex;
                Scenario(row) = cfg.scenarioNames(scenarioIndex);
                Separation(row) = cfg.separationValues(dIndex);
                Metabolite(row) = metaboliteIndex;
                Valid(row) = sum(codes == 1);
                NonpositiveDenominator(row) = sum(codes == 2);
                NonfiniteQ(row) = sum(codes == 3);
                QOutsideOpenUnitInterval(row) = sum(codes == 4);
            end
        end
    end
    status = table(ScenarioIndex, Scenario, Separation, Metabolite, ...
        Valid, NonpositiveDenominator, NonfiniteQ, ...
        QOutsideOpenUnitInterval);
end

function [bias, rmse, medianAbsoluteError] = errorMetrics(estimate, truth)
    if isempty(estimate)
        bias = NaN;
        rmse = NaN;
        medianAbsoluteError = NaN;
        return;
    end
    error = estimate - truth;
    bias = mean(error);
    rmse = sqrt(mean(error.^2));
    medianAbsoluteError = median(abs(error));
end

function value = meanForMask(array, mask)
    array = reshape(array, [], 1);
    mask = reshape(mask, [], 1) & isfinite(array);
    if any(mask)
        value = mean(array(mask));
    else
        value = NaN;
    end
end
