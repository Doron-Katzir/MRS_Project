%% TwoGaussianOverlapCovarianceSimulation
% Compare covariance estimators in a two-Gaussian repeated-MRS toy model.
%
% For repetition t, the simulated spectrum and concentration vector obey
%
%   s_t = B*c_t + epsilon_t,
%   c_t ~ N(mu_true, Sigma_c_true),
%   epsilon_t ~ N(0, sigma2_true*I_K).
%
% The Gaussian peak separation is varied through
%
%   d = abs(center2 - center1)/peakWidth.
%
% Small d means that the basis functions are difficult to separate. This
% overlap sweep tests whether joint marginal maximum likelihood estimates
% the true concentration covariance better as the basis functions become
% less separable.
%
% Three concentration-covariance estimators are compared:
%   1. Naive OLS covariance.
%   2. Noise-corrected OLS covariance.
%   3. Joint ordinary marginal maximum likelihood (ML), not REML.
%
% Naive OLS estimates Sigma_c + sigma2*(B'*B)^(-1), rather than Sigma_c
% alone. Corrected OLS subtracts an estimate of this fitting-noise
% covariance. Joint ML instead fits the spectra directly through
% V = B*Sigma_c*B' + sigma2*I_K.

clear;
clc;
close all;

%% Simulation settings

K = 160;
T = 36;
x = linspace(0, 10, K).';
peakWidth = 0.55;
dValues = [0.25, 0.5, 0.75, 1, 1.5, 2, 3, 4].';
nMC = 500;
randomSeed = 1;

% Optional environment overrides support small automated smoke tests while
% leaving the scientific defaults above unchanged for ordinary use.
nMCOverride = str2double(getenv("TWO_GAUSSIAN_NMC"));
if isfinite(nMCOverride) && nMCOverride >= 1 ...
        && nMCOverride == floor(nMCOverride)
    nMC = nMCOverride;
end
dValuesOverrideText = strrep(getenv("TWO_GAUSSIAN_DVALUES"), ",", " ");
dValuesOverride = sscanf(dValuesOverrideText, "%f").';
if ~isempty(dValuesOverride) && all(isfinite(dValuesOverride)) ...
        && all(dValuesOverride > 0)
    dValues = dValuesOverride(:);
end

mu_true = [3.0; 2.0];
tau1 = 0.65;
tau2 = 0.45;
rho_c_true = 0.4;
tau12 = rho_c_true * tau1 * tau2;
Sigma_c_true = [tau1^2, tau12; tau12, tau2^2];
sigma_true = 0.9;
sigma2_true = sigma_true^2;

rng(randomSeed, "twister");

% Write outputs next to this script, independent of MATLAB's current folder.
scriptPath = mfilename("fullpath");
if strlength(scriptPath) == 0
    outputDir = pwd;
else
    outputDir = fileparts(scriptPath);
end
outputDirOverride = string(getenv("TWO_GAUSSIAN_OUTPUT_DIR"));
if strlength(outputDirOverride) > 0
    outputDir = outputDirOverride;
end
if ~isfolder(outputDir)
    mkdir(outputDir);
end

methodNames = ["NaiveOLS", "CorrectedOLS", "JointML"];
nMethods = numel(methodNames);
nD = numel(dValues);

%% Preallocate the full Monte Carlo result structure

rawResults = struct();
rawResults.methodNames = methodNames;
rawResults.dValues = dValues;
rawResults.basisOverlapRho = nan(nD, 1);
rawResults.condGram = nan(nD, 1);

% Metric arrays are nD-by-nMC-by-nMethods.
rawResults.frobeniusError = nan(nD, nMC, nMethods);
rawResults.tau1sqError = nan(nD, nMC, nMethods);
rawResults.tau2sqError = nan(nD, nMC, nMethods);
rawResults.tau12Error = nan(nD, nMC, nMethods);
rawResults.rhoHat = nan(nD, nMC, nMethods);
rawResults.rhoError = nan(nD, nMC, nMethods);

% Store the underlying estimates too, not only their derived metrics.
% Covariance arrays are 2-by-2-by-nMC-by-nD.
rawResults.SigmaNaiveOLS = nan(2, 2, nMC, nD);
rawResults.SigmaCorrectedOLS = nan(2, 2, nMC, nD);
rawResults.SigmaJointML = nan(2, 2, nMC, nD);
rawResults.muOLS = nan(2, nMC, nD);
rawResults.muML = nan(2, nMC, nD);

rawResults.sigma2OLS = nan(nD, nMC);
rawResults.sigma2ML = nan(nD, nMC);
rawResults.correctedOLSIsPSD = false(nD, nMC);
rawResults.mlSuccess = false(nD, nMC);
rawResults.mlExitFlag = nan(nD, nMC);
rawResults.mlFinalNLL = nan(nD, nMC);
rawResults.mlOptimizer = strings(nD, nMC);

%% Monte Carlo simulation across basis overlaps

for dIndex = 1:nD
    d = dValues(dIndex);
    [B, b1, b2] = makeGaussianBasis(x, d, peakWidth);

    gramMatrix = B.' * B;
    rawResults.basisOverlapRho(dIndex) = (b1.' * b2) / ...
        sqrt((b1.' * b1) * (b2.' * b2));
    rawResults.condGram(dIndex) = cond(gramMatrix);

    fprintf("d = %.2f (%d/%d): basis overlap rho = %.6f, " + ...
        "cond(B''B) = %.4g\n", d, dIndex, nD, ...
        rawResults.basisOverlapRho(dIndex), ...
        rawResults.condGram(dIndex));

    for mcIndex = 1:nMC
        [S, ~] = simulateDataset(B, T, mu_true, ...
            Sigma_c_true, sigma2_true);

        olsFit = fitOLS(S, B);
        rawResults.SigmaNaiveOLS(:, :, mcIndex, dIndex) = ...
            olsFit.SigmaNaive;
        rawResults.SigmaCorrectedOLS(:, :, mcIndex, dIndex) = ...
            olsFit.SigmaCorrected;
        rawResults.muOLS(:, mcIndex, dIndex) = olsFit.mu;
        rawResults.sigma2OLS(dIndex, mcIndex) = olsFit.sigma2;
        rawResults.correctedOLSIsPSD(dIndex, mcIndex) = ...
            olsFit.correctedIsPSD;

        naiveMetrics = computeMetrics( ...
            olsFit.SigmaNaive, Sigma_c_true, rho_c_true);
        correctedMetrics = computeMetrics( ...
            olsFit.SigmaCorrected, Sigma_c_true, rho_c_true);
        rawResults = storeMetrics(rawResults, naiveMetrics, ...
            dIndex, mcIndex, 1);
        rawResults = storeMetrics(rawResults, correctedMetrics, ...
            dIndex, mcIndex, 2);

        % Any numerical failure is recorded and the simulation continues.
        try
            mlFit = fitJointML(S, B, olsFit);
        catch fitException
            mlFit = failedMLFit(fitException.message);
        end

        rawResults.SigmaJointML(:, :, mcIndex, dIndex) = mlFit.Sigma;
        rawResults.muML(:, mcIndex, dIndex) = mlFit.mu;
        rawResults.sigma2ML(dIndex, mcIndex) = mlFit.sigma2;
        rawResults.mlSuccess(dIndex, mcIndex) = mlFit.success;
        rawResults.mlExitFlag(dIndex, mcIndex) = mlFit.exitFlag;
        rawResults.mlFinalNLL(dIndex, mcIndex) = mlFit.finalNLL;
        rawResults.mlOptimizer(dIndex, mcIndex) = mlFit.optimizer;

        % Exclude unsuccessful optimizations from ML accuracy summaries.
        % Their frequency is retained in the reported ML failure rate.
        if mlFit.success
            mlMetrics = computeMetrics( ...
                mlFit.Sigma, Sigma_c_true, rho_c_true);
            rawResults = storeMetrics(rawResults, mlMetrics, ...
                dIndex, mcIndex, 3);
        end
    end

    fprintf("  Completed %d Monte Carlo repetitions; ML failures: %d.\n", ...
        nMC, sum(~rawResults.mlSuccess(dIndex, :)));
end

%% Summarize, export, and save full results

summaryTable = summarizeResults(dValues, ...
    rawResults.basisOverlapRho, rawResults.condGram, rawResults);

summaryCsvPath = fullfile(outputDir, "TwoGaussianOverlap_Summary.csv");
rawMatPath = fullfile(outputDir, "TwoGaussianOverlap_RawResults.mat");
writetable(summaryTable, summaryCsvPath);

simulationConfig = struct( ...
    "K", K, ...
    "T", T, ...
    "x", x, ...
    "peakWidth", peakWidth, ...
    "dValues", dValues, ...
    "nMC", nMC, ...
    "randomSeed", randomSeed, ...
    "mu_true", mu_true, ...
    "Sigma_c_true", Sigma_c_true, ...
    "rho_c_true", rho_c_true, ...
    "sigma_true", sigma_true, ...
    "sigma2_true", sigma2_true);

save(rawMatPath, "rawResults", "summaryTable", "simulationConfig", "-v7.3");

%% Figures

lineWidth = 1.8;
markerSize = 7;
plotColors = lines(nMethods);

fig = figure("Name", "Frobenius error versus overlap", "Color", "w");
hold on;
plot(dValues, summaryTable.MeanFrobeniusError_NaiveOLS, "-o", ...
    "Color", plotColors(1, :), "LineWidth", lineWidth, ...
    "MarkerSize", markerSize, "DisplayName", "Naive OLS");
plot(dValues, summaryTable.MeanFrobeniusError_CorrectedOLS, "-s", ...
    "Color", plotColors(2, :), "LineWidth", lineWidth, ...
    "MarkerSize", markerSize, "DisplayName", "Corrected OLS");
plot(dValues, summaryTable.MeanFrobeniusError_JointML, "-^", ...
    "Color", plotColors(3, :), "LineWidth", lineWidth, ...
    "MarkerSize", markerSize, "DisplayName", "Joint ML");
xlabel("Normalized peak separation d");
ylabel("Mean Frobenius error");
title("Concentration-covariance error versus Gaussian separation");
legend("Location", "best");
grid on;
saveFigurePair(fig, outputDir, "Fig_FrobeniusError_vs_Overlap");

fig = figure("Name", "Tau12 RMSE versus overlap", "Color", "w");
hold on;
plot(dValues, summaryTable.RMSE_Tau12_NaiveOLS, "-o", ...
    "Color", plotColors(1, :), "LineWidth", lineWidth, ...
    "MarkerSize", markerSize, "DisplayName", "Naive OLS");
plot(dValues, summaryTable.RMSE_Tau12_CorrectedOLS, "-s", ...
    "Color", plotColors(2, :), "LineWidth", lineWidth, ...
    "MarkerSize", markerSize, "DisplayName", "Corrected OLS");
plot(dValues, summaryTable.RMSE_Tau12_JointML, "-^", ...
    "Color", plotColors(3, :), "LineWidth", lineWidth, ...
    "MarkerSize", markerSize, "DisplayName", "Joint ML");
xlabel("Normalized peak separation d");
ylabel("RMSE of covariance \tau_{12}");
title("Cross-covariance estimation error versus Gaussian separation");
legend("Location", "best");
grid on;
saveFigurePair(fig, outputDir, "Fig_Tau12RMSE_vs_Overlap");

fig = figure("Name", "Concentration-correlation RMSE versus overlap", ...
    "Color", "w");
hold on;
plot(dValues, summaryTable.RMSE_RhoC_NaiveOLS, "-o", ...
    "Color", plotColors(1, :), "LineWidth", lineWidth, ...
    "MarkerSize", markerSize, "DisplayName", "Naive OLS");
plot(dValues, summaryTable.RMSE_RhoC_CorrectedOLS, "-s", ...
    "Color", plotColors(2, :), "LineWidth", lineWidth, ...
    "MarkerSize", markerSize, "DisplayName", "Corrected OLS");
plot(dValues, summaryTable.RMSE_RhoC_JointML, "-^", ...
    "Color", plotColors(3, :), "LineWidth", lineWidth, ...
    "MarkerSize", markerSize, "DisplayName", "Joint ML");
xlabel("Normalized peak separation d");
ylabel("RMSE of concentration correlation \rho_c");
title("Concentration-correlation error versus Gaussian separation");
legend("Location", "best");
grid on;
saveFigurePair(fig, outputDir, "Fig_RhoRMSE_vs_Overlap");

fig = figure("Name", "Spectral-noise variance versus overlap", ...
    "Color", "w");
hold on;
yline(sigma2_true, "w--", "LineWidth", 1.5, ...
    "DisplayName", "True \sigma^2");
plot(dValues, summaryTable.MeanSigma2_OLS, "-o", ...
    "Color", plotColors(1, :), "LineWidth", lineWidth, ...
    "MarkerSize", markerSize, "DisplayName", "OLS residual estimate");
plot(dValues, summaryTable.MeanSigma2_ML, "-^", ...
    "Color", plotColors(3, :), "LineWidth", lineWidth, ...
    "MarkerSize", markerSize, "DisplayName", "Joint ML");
xlabel("Normalized peak separation d");
ylabel("Mean estimate of \sigma^2");
title("Spectral-noise variance estimation");
legend("Location", "best");
grid on;
saveFigurePair(fig, outputDir, "Fig_Sigma2_vs_Overlap");

fig = figure("Name", "Corrected OLS non-PSD rate", "Color", "w");
plot(dValues, summaryTable.CorrectedOLS_NonPSDRate, "-o", ...
    "Color", plotColors(2, :), "LineWidth", lineWidth, ...
    "MarkerSize", markerSize);
xlabel("Normalized peak separation d");
ylabel("Non-PSD rate");
ylim([0, 1]);
title("Invalid corrected-OLS covariance estimates");
grid on;
saveFigurePair(fig, outputDir, "Fig_OLSNonPSD_vs_Overlap");

fprintf("\nMain summary table:\n");
disp(summaryTable);
fprintf("Summary CSV: %s\n", summaryCsvPath);
fprintf("Raw MAT file: %s\n", rawMatPath);

%% Local functions

function [B, b1, b2] = makeGaussianBasis(x, d, peakWidth)
%MAKEGAUSSIANBASIS Construct two max-normalized Gaussian basis vectors.

    center1 = 5 - (d * peakWidth) / 2;
    center2 = 5 + (d * peakWidth) / 2;

    b1 = exp(-(x - center1).^2 / (2 * peakWidth^2));
    b2 = exp(-(x - center2).^2 / (2 * peakWidth^2));
    b1 = b1 / max(b1);
    b2 = b2 / max(b2);
    B = [b1, b2];
end

function [S, C_true] = simulateDataset(B, T, mu_true, ...
        Sigma_c_true, sigma2_true)
%SIMULATEDATASET Draw concentrations and spectra from the stated model.
%   This construction avoids requiring the Statistics and Machine Learning
%   Toolbox for multivariate-normal random draws.

    concentrationChol = chol(Sigma_c_true, "lower");
    C_true = repmat(mu_true.', T, 1) ...
        + randn(T, 2) * concentrationChol.';
    E = sqrt(sigma2_true) * randn(T, size(B, 1));
    S = C_true * B.' + E;
end

function fit = fitOLS(S, B)
%FITOLS Fit each spectrum separately and form both OLS covariances.

    [T, K] = size(S);
    nBasis = size(B, 2);
    gramMatrix = B.' * B;

    C_OLS = (gramMatrix \ (B.' * S.')).';
    mu_OLS = mean(C_OLS, 1).';
    Sigma_OLS_naive = cov(C_OLS, 1); % ML divisor T, not T-1.
    Sigma_OLS_naive = (Sigma_OLS_naive + Sigma_OLS_naive.') / 2;

    S_fit_OLS = C_OLS * B.';
    R_OLS = S - S_fit_OLS;
    sigma2_OLS = sum(R_OLS(:).^2) / (T * (K - nBasis));

    % Use a solve rather than explicitly forming inv(B'*B).
    Sigma_OLS_corr = Sigma_OLS_naive ...
        - sigma2_OLS * (gramMatrix \ eye(nBasis));
    Sigma_OLS_corr = (Sigma_OLS_corr + Sigma_OLS_corr.') / 2;
    correctedEigenvalues = eig(Sigma_OLS_corr);

    fit.C = C_OLS;
    fit.mu = mu_OLS;
    fit.SigmaNaive = Sigma_OLS_naive;
    fit.SigmaCorrected = Sigma_OLS_corr;
    fit.sigma2 = sigma2_OLS;
    fit.correctedIsPSD = all(correctedEigenvalues >= -1e-10);
end

function fit = fitJointML(S, B, olsFit)
%FITJOINTML Fit mu, Sigma_c, and sigma2 by ordinary marginal ML.
%   theta = [mu1, mu2, a, b, r, dNoise], with
%       L = [exp(a), 0; r, exp(b)], Sigma_c = L*L',
%       sigma2 = exp(dNoise).

    sigma2_0 = max(olsFit.sigma2, 1e-10);
    Sigma0_PD = makePositiveDefiniteForInitialization( ...
        olsFit.SigmaCorrected, olsFit.SigmaNaive);
    L0 = chol(Sigma0_PD, "lower");

    theta0 = [olsFit.mu(1), olsFit.mu(2), ...
        log(L0(1, 1)), log(L0(2, 2)), L0(2, 1), log(sigma2_0)];
    objective = @(theta) marginalNegativeLogLikelihood(theta, S, B);

    optimizer = "fminsearch";
    if exist("fminunc", "file") == 2
        optimizer = "fminunc";
        options = optimoptions("fminunc", ...
            "Display", "off", ...
            "Algorithm", "quasi-newton", ...
            "MaxIterations", 500, ...
            "MaxFunctionEvaluations", 3000);
        try
            [thetaHat, finalNLL, exitFlag] = ...
                fminunc(objective, theta0, options);
        catch
            % A file may exist even when its toolbox license is unavailable.
            optimizer = "fminsearch-fallback";
            options = optimset( ...
                "Display", "off", ...
                "MaxIter", 1000, ...
                "MaxFunEvals", 5000, ...
                "TolX", 1e-8, ...
                "TolFun", 1e-8);
            [thetaHat, finalNLL, exitFlag] = ...
                fminsearch(objective, theta0, options);
        end
    else
        options = optimset( ...
            "Display", "off", ...
            "MaxIter", 1000, ...
            "MaxFunEvals", 5000, ...
            "TolX", 1e-8, ...
            "TolFun", 1e-8);
        [thetaHat, finalNLL, exitFlag] = ...
            fminsearch(objective, theta0, options);
    end

    [muHat, SigmaHat, sigma2Hat] = unpackTheta(thetaHat);
    SigmaHat = (SigmaHat + SigmaHat.') / 2;
    finiteFit = all(isfinite(thetaHat)) && isfinite(finalNLL) ...
        && all(isfinite(SigmaHat), "all") && isfinite(sigma2Hat);

    fit.mu = muHat;
    fit.Sigma = SigmaHat;
    fit.sigma2 = sigma2Hat;
    fit.theta = thetaHat;
    fit.exitFlag = exitFlag;
    fit.finalNLL = finalNLL;
    fit.optimizer = optimizer;
    fit.success = exitFlag > 0 && finiteFit && sigma2Hat > 0;
end

function nll = marginalNegativeLogLikelihood(theta, S, B)
%MARGINALNEGATIVELOGLIKELIHOOD Evaluate the ordinary marginal Gaussian NLL.
%   Constants independent of theta are omitted. The Cholesky factor of the
%   full spectral covariance is used for both the log determinant and the
%   quadratic form; inv(V) is never formed.

    largePenalty = 1e100;
    if numel(theta) ~= 6 || any(~isfinite(theta))
        nll = largePenalty;
        return;
    end

    [mu, Sigma_c, sigma2] = unpackTheta(theta);
    if any(~isfinite(Sigma_c), "all") || ~isfinite(sigma2) || sigma2 <= 0
        nll = largePenalty;
        return;
    end

    K = size(B, 1);
    T = size(S, 1);
    V = B * Sigma_c * B.' + sigma2 * eye(K);
    V = (V + V.') / 2;
    [Lchol, cholFailure] = chol(V, "lower");
    if cholFailure ~= 0 || any(~isfinite(Lchol), "all")
        nll = largePenalty;
        return;
    end

    residuals = S.' - B * mu; % K-by-T.
    alpha = Lchol \ residuals;
    quad = sum(alpha(:).^2);
    logdetV = 2 * sum(log(diag(Lchol)));
    nll = 0.5 * T * logdetV + 0.5 * quad;

    if ~isfinite(nll)
        nll = largePenalty;
    end
end

function [mu, Sigma_c, sigma2, L] = unpackTheta(theta)
%UNPACKTHETA Convert unconstrained optimizer parameters to model parameters.

    mu = theta(1:2).';
    L = [exp(theta(3)), 0; theta(5), exp(theta(4))];
    Sigma_c = L * L.';
    sigma2 = exp(theta(6));
end

function SigmaPD = makePositiveDefiniteForInitialization( ...
        Sigma0, SigmaNaive)
%MAKEPOSITIVEDEFINITEFORINITIALIZATION Repair only the ML starting value.
%   The reported corrected-OLS estimate itself is never altered.

    Sigma0 = (Sigma0 + Sigma0.') / 2;
    SigmaNaive = (SigmaNaive + SigmaNaive.') / 2;
    varianceScale = max(trace(SigmaNaive) / 2, eps);
    eigFloor = max(1e-6 * varianceScale, 1e-10);

    [eigenvectors, eigenvalueMatrix] = eig(Sigma0);
    eigenvalues = diag(eigenvalueMatrix);
    eigenvalues(~isfinite(eigenvalues)) = eigFloor;
    eigenvalues = max(eigenvalues, eigFloor);
    SigmaPD = eigenvectors * diag(eigenvalues) * eigenvectors.';
    SigmaPD = (SigmaPD + SigmaPD.') / 2;

    [~, cholFailure] = chol(SigmaPD, "lower");
    if cholFailure ~= 0 || any(~isfinite(SigmaPD), "all")
        fallbackDiagonal = max(diag(SigmaNaive), eigFloor);
        fallbackDiagonal(~isfinite(fallbackDiagonal)) = eigFloor;
        SigmaPD = diag(fallbackDiagonal);
    end
end

function metrics = computeMetrics(SigmaHat, SigmaTrue, rhoTrue)
%COMPUTEMETRICS Calculate covariance and concentration-correlation errors.

    SigmaHat = (SigmaHat + SigmaHat.') / 2;
    metrics.frobeniusError = norm(SigmaHat - SigmaTrue, "fro");
    metrics.tau1sqError = SigmaHat(1, 1) - SigmaTrue(1, 1);
    metrics.tau2sqError = SigmaHat(2, 2) - SigmaTrue(2, 2);
    metrics.tau12Error = SigmaHat(1, 2) - SigmaTrue(1, 2);

    denominatorSquared = SigmaHat(1, 1) * SigmaHat(2, 2);
    if denominatorSquared > 0 && isfinite(denominatorSquared)
        metrics.rhoHat = SigmaHat(1, 2) / sqrt(denominatorSquared);
        metrics.rhoError = metrics.rhoHat - rhoTrue;
    else
        metrics.rhoHat = NaN;
        metrics.rhoError = NaN;
    end
end

function summaryTable = summarizeResults(dValues, basisOverlapRho, ...
        condGram, rawResults)
%SUMMARIZERESULTS Aggregate Monte Carlo errors at each overlap value.

    nD = numel(dValues);
    nMethods = size(rawResults.frobeniusError, 3);

    meanFrobenius = nan(nD, nMethods);
    rmseTau1Sq = nan(nD, nMethods);
    rmseTau2Sq = nan(nD, nMethods);
    rmseTau12 = nan(nD, nMethods);
    rmseRhoC = nan(nD, nMethods);

    for dIndex = 1:nD
        for methodIndex = 1:nMethods
            meanFrobenius(dIndex, methodIndex) = finiteMean( ...
                rawResults.frobeniusError(dIndex, :, methodIndex));
            rmseTau1Sq(dIndex, methodIndex) = finiteRMSE( ...
                rawResults.tau1sqError(dIndex, :, methodIndex));
            rmseTau2Sq(dIndex, methodIndex) = finiteRMSE( ...
                rawResults.tau2sqError(dIndex, :, methodIndex));
            rmseTau12(dIndex, methodIndex) = finiteRMSE( ...
                rawResults.tau12Error(dIndex, :, methodIndex));
            rmseRhoC(dIndex, methodIndex) = finiteRMSE( ...
                rawResults.rhoError(dIndex, :, methodIndex));
        end
    end

    meanSigma2OLS = nan(nD, 1);
    meanSigma2ML = nan(nD, 1);
    for dIndex = 1:nD
        meanSigma2OLS(dIndex) = finiteMean( ...
            rawResults.sigma2OLS(dIndex, :));
        successfulSigma2 = rawResults.sigma2ML(dIndex, ...
            rawResults.mlSuccess(dIndex, :));
        meanSigma2ML(dIndex) = finiteMean(successfulSigma2);
    end

    correctedOLSNonPSDRate = mean(~rawResults.correctedOLSIsPSD, 2);
    mlFailureRate = mean(~rawResults.mlSuccess, 2);

    summaryTable = table( ...
        dValues, basisOverlapRho, condGram, ...
        meanFrobenius(:, 1), meanFrobenius(:, 2), meanFrobenius(:, 3), ...
        rmseTau1Sq(:, 1), rmseTau1Sq(:, 2), rmseTau1Sq(:, 3), ...
        rmseTau2Sq(:, 1), rmseTau2Sq(:, 2), rmseTau2Sq(:, 3), ...
        rmseTau12(:, 1), rmseTau12(:, 2), rmseTau12(:, 3), ...
        rmseRhoC(:, 1), rmseRhoC(:, 2), rmseRhoC(:, 3), ...
        meanSigma2OLS, meanSigma2ML, correctedOLSNonPSDRate, ...
        mlFailureRate);
    summaryTable.Properties.VariableNames = [ ...
        "d", "BasisOverlapRho", "CondGram", ...
        "MeanFrobeniusError_NaiveOLS", ...
        "MeanFrobeniusError_CorrectedOLS", ...
        "MeanFrobeniusError_JointML", ...
        "RMSE_Tau1Sq_NaiveOLS", "RMSE_Tau1Sq_CorrectedOLS", ...
        "RMSE_Tau1Sq_JointML", ...
        "RMSE_Tau2Sq_NaiveOLS", "RMSE_Tau2Sq_CorrectedOLS", ...
        "RMSE_Tau2Sq_JointML", ...
        "RMSE_Tau12_NaiveOLS", "RMSE_Tau12_CorrectedOLS", ...
        "RMSE_Tau12_JointML", ...
        "RMSE_RhoC_NaiveOLS", "RMSE_RhoC_CorrectedOLS", ...
        "RMSE_RhoC_JointML", ...
        "MeanSigma2_OLS", "MeanSigma2_ML", ...
        "CorrectedOLS_NonPSDRate", "ML_FailureRate"];
end

function rawResults = storeMetrics(rawResults, metrics, ...
        dIndex, mcIndex, methodIndex)
%STOREMETRICS Place one method's scalar metrics in preallocated arrays.

    rawResults.frobeniusError(dIndex, mcIndex, methodIndex) = ...
        metrics.frobeniusError;
    rawResults.tau1sqError(dIndex, mcIndex, methodIndex) = ...
        metrics.tau1sqError;
    rawResults.tau2sqError(dIndex, mcIndex, methodIndex) = ...
        metrics.tau2sqError;
    rawResults.tau12Error(dIndex, mcIndex, methodIndex) = ...
        metrics.tau12Error;
    rawResults.rhoHat(dIndex, mcIndex, methodIndex) = metrics.rhoHat;
    rawResults.rhoError(dIndex, mcIndex, methodIndex) = metrics.rhoError;
end

function fit = failedMLFit(~)
%FAILEDMLFIT Return a standardized failed-fit record after an exception.

    fit.mu = nan(2, 1);
    fit.Sigma = nan(2, 2);
    fit.sigma2 = NaN;
    fit.theta = nan(1, 6);
    fit.exitFlag = -999;
    fit.finalNLL = NaN;
    fit.optimizer = "error";
    fit.success = false;
end

function value = finiteMean(values)
%FINITEMEAN Mean after excluding NaN and infinite values.

    values = values(isfinite(values));
    if isempty(values)
        value = NaN;
    else
        value = mean(values);
    end
end

function value = finiteRMSE(errors)
%FINITERMSE Root-mean-square error after excluding invalid values.

    errors = errors(isfinite(errors));
    if isempty(errors)
        value = NaN;
    else
        value = sqrt(mean(errors.^2));
    end
end

function saveFigurePair(fig, outputDir, baseName)
%SAVEFIGUREPAIR Save one figure in both PNG and editable FIG formats.

    exportgraphics(fig, fullfile(outputDir, baseName + ".png"), ...
        "Resolution", 200);
    savefig(fig, fullfile(outputDir, baseName + ".fig"));
end
