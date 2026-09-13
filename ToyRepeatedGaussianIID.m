%% ToyRepeatedGaussianIID
% A transparent simulation and fit for repeated spectra with one known basis:
%
%   s_t = c_t b + epsilon_t,
%   c_t ~ N(mu, tau^2),       epsilon_t ~ N(0, sigma^2 I_K).
%
% After integrating out c_t, the repeated spectra are IID with
%
%   s_t ~ N(mu*b, tau^2*b*b' + sigma^2*I_K).
%
% This script estimates mu, tau^2, and sigma^2 jointly, then compares the
% independent OLS concentration estimates with the fitted-model posterior
% means. All spectra are real-valued in this first toy problem.

clear;
clc;
close all;

%% 1. Simulation settings

K = 160;                  % Number of samples in each spectrum
T = 36;                   % Number of repeated measurements
xMin = 0;
xMax = 10;
gaussianCenter = 5.0;
gaussianWidth = 0.55;     % Standard-deviation width in x-axis units

muTrue = 3.0;
tauTrue = 0.65;           % Population standard deviation of concentration
sigmaTrue = 0.9;         % Noise standard deviation at each spectral sample

randomSeed = 7;
doLikelihoodCheck = true; % Optional one-dimensional likelihood diagnostic

% Monte Carlo settings used in Section 12.
nMonteCarlo = 1000;
monteCarloSeed = 123;

rng(randomSeed, "twister");

%% 2. Generate the known Gaussian basis, b (K x 1)

x = linspace(xMin, xMax, K).';
b = exp(-(x - gaussianCenter).^2 / (2 * gaussianWidth^2));

% Normalize max(b) to one. Consequently, c_t is the noiseless peak height.
b = b / max(b);

figure("Name", "Known Gaussian basis", "Color", "w");
plot(x, b, "LineWidth", 2);
xlabel("Spectral coordinate x");
ylabel("b(x)");
title("Known Gaussian basis, normalized so max(b) = 1");
grid on;

%% 3. Generate the true concentrations, cTrue (T x 1)

cTrue = muTrue + tauTrue * randn(T, 1);
measurementIndex = (1:T).';

figure("Name", "True concentrations", "Color", "w");
plot(measurementIndex, cTrue, "o-", "LineWidth", 1.2, ...
    "MarkerSize", 4);
yline(muTrue, "--", "True population mean", "LineWidth", 1.5);
xlabel("Measurement index t");
ylabel("Concentration c_t");
title("Simulated IID concentrations");
grid on;

%% 4. Generate repeated noisy spectra

% S is K x T: column t is the measured spectrum s_t.
noise = sigmaTrue * randn(K, T);
S = b * cTrue.' + noise;

figure("Name", "Repeated noisy spectra", "Color", "w");
hSpectra = plot(x, S, "Color", [0.72, 0.80, 0.92], "LineWidth", 0.6);
hold on;
hMean = plot(x, mean(S, 2), "k", "LineWidth", 2.2);
xlabel("Spectral coordinate x");
ylabel("Measured signal");
title(sprintf("All %d repeated spectra", T));
legend([hSpectra(1), hMean], "Individual spectra", "Measured mean", ...
    "Location", "best");
grid on;

figure("Name", "Mean measured spectrum", "Color", "w");
plot(x, mean(S, 2), "k", "LineWidth", 2.2, "Color", "w");
hold on;
plot(x, muTrue * b, "r--", "LineWidth", 2);
xlabel("Spectral coordinate x");
ylabel("Signal");
title("Mean measured spectrum and true population mean");
legend("Mean measured spectrum", "\mu_{true} b", "Location", "best");
grid on;

%% 5-6. Fit the marginal Gaussian likelihood

% The optimizer uses theta = [mu, log(tau^2), log(sigma^2)]. Taking the
% exponential inside the objective ensures that both variances are positive.
%
% Simple data-based starting values improve fminsearch convergence without
% hiding any model assumptions. FitMarginalGaussianIID constructs these
% values and is also used unchanged by the Monte Carlo experiment below.
bTb = b.' * b;

optimizerOptions = optimset( ...
    "Display", "final", ...
    "MaxFunEvals", 3000, ...
    "MaxIter", 1500, ...
    "TolX", 1e-9, ...
    "TolFun", 1e-9);

fitSingle = FitMarginalGaussianIID(S, b, optimizerOptions);

thetaHat = fitSingle.thetaHat;
nllHat = fitSingle.negativeLogLikelihood;
exitFlag = fitSingle.exitFlag;
optimizerOutput = fitSingle.optimizerOutput;
cOLS = fitSingle.cOLS;                    % Also used for comparison below
thetaInitial = fitSingle.thetaInitial;
muInitial = thetaInitial(1);
tau2Initial = exp(thetaInitial(2));
sigma2Initial = exp(thetaInitial(3));

% Reuse the exact same objective for the optional likelihood plot.
objective = @(theta) MarginalNegativeLogLikelihood( ...
    theta, fitSingle.likelihoodStatistics);

muHat = thetaHat(1);
tau2Hat = exp(thetaHat(2));
sigma2Hat = exp(thetaHat(3));
tauHat = sqrt(tau2Hat);
sigmaHat = sqrt(sigma2Hat);

parameter = ["mu"; "tau^2"; "tau"; "sigma^2"; "sigma"];
trueValue = [muTrue; tauTrue^2; tauTrue; sigmaTrue^2; sigmaTrue];
estimatedValue = [muHat; tau2Hat; tauHat; sigma2Hat; sigmaHat];
absoluteError = abs(estimatedValue - trueValue);
parameterComparison = table(parameter, trueValue, estimatedValue, absoluteError);

fprintf("\nPopulation parameter comparison:\n");
disp(parameterComparison);
fprintf("Final negative log-likelihood = %.6f\n", nllHat);
fprintf("fminsearch exit flag = %d (%d iterations)\n", ...
    exitFlag, optimizerOutput.iterations);

%% 7. Independent per-spectrum estimates

% Ordinary least squares treats every spectrum separately:
%
%   c_t^OLS = (b' s_t) / (b' b).
%
% cOLS was calculated above and is stored as a T x 1 vector.

%% 8. Individual posterior means under the fitted IID model

% Gaussian conditioning gives
%
% E[c_t | s_t] = muHat ...
%   + tau2Hat*b'*(s_t - muHat*b)/(sigma2Hat + tau2Hat*b'*b).
%
% The multiplier shrinks the noisy OLS information toward the fitted
% population mean. cJoint is T x 1.
centeredSpectra = S - muHat * b;
posteriorGain = tau2Hat / (sigma2Hat + tau2Hat * bTb);
cJoint = (muHat + posteriorGain * (b.' * centeredSpectra)).';

% The posterior variance is identical for all t in this IID model.
posteriorVariance = 1 / (1 / tau2Hat + bTb / sigma2Hat);

%% 9. Compare individual concentration estimates

figure("Name", "Concentration estimates", "Color", "w");
plot(measurementIndex, cTrue, "ko-", "LineWidth", 1.4, ...
    "MarkerSize", 4, "DisplayName", "True c_t", "Color", "y");
hold on;
plot(measurementIndex, cOLS, ".-", "Color", [0.85, 0.33, 0.10], ...
    "LineWidth", 1, "MarkerSize", 12, "DisplayName", "Independent OLS");
plot(measurementIndex, cJoint, ".-", "Color", [0, 0.45, 0.74], ...
    "LineWidth", 1, "MarkerSize", 12, "DisplayName", "IID posterior mean");
xlabel("Measurement index t");
ylabel("Concentration");
title("True and estimated concentrations");
legend("Location", "best");
grid on;

figure("Name", "True versus estimated concentration", "Color", "w");
scatter(cTrue, cOLS, 42, [0.85, 0.33, 0.10], "o", ...
    "DisplayName", "Independent OLS");
hold on;
scatter(cTrue, cJoint, 42, [0, 0.45, 0.74], "filled", ...
    "DisplayName", "IID posterior mean");
identityLimits = [min([cTrue; cOLS; cJoint]), max([cTrue; cOLS; cJoint])];
plot(identityLimits, identityLimits, "k--", "LineWidth", 1.5, ...
    "DisplayName", "Identity line");
xlim(identityLimits);
ylim(identityLimits);
axis square;
xlabel("True concentration");
ylabel("Estimated concentration");
title("Concentration recovery");
legend("Location", "best");
grid on;

%% 10. Numerical performance summary

rmseOLS = sqrt(mean((cTrue - cOLS).^2));
rmseJoint = sqrt(mean((cTrue - cJoint).^2));

fprintf("\nConcise simulation summary\n");
fprintf("True mu          = %.6f\n", muTrue);
fprintf("Estimated mu     = %.6f\n\n", muHat);
fprintf("True tau         = %.6f\n", tauTrue);
fprintf("Measured tau     = %.6f\n", std(cTrue,1));
fprintf("Estimated tau    = %.6f\n\n", tauHat);
fprintf("True sigma       = %.6f\n", sigmaTrue);
fprintf("Estimated sigma  = %.6f\n\n", sigmaHat);
fprintf("OLS RMSE         = %.6f\n", rmseOLS);
fprintf("Joint RMSE       = %.6f\n", rmseJoint);
fprintf("Posterior SD     = %.6f (same for every measurement)\n", ...
    sqrt(posteriorVariance));

%% 11. Optional likelihood sanity check

if doLikelihoodCheck
    % Hold the fitted variances fixed and vary mu near its optimum.
    muHalfWidth = max(3 * tauHat / sqrt(T), 0.25);
    muGrid = linspace(muHat - muHalfWidth, muHat + muHalfWidth, 121);
    nllAlongMu = arrayfun(@(mu) objective( ...
        [mu, log(tau2Hat), log(sigma2Hat)]), muGrid);

    figure("Name", "Likelihood sanity check", "Color", "w");
    plot(muGrid, nllAlongMu - nllHat, "LineWidth", 2);
    hold on;
    xline(muHat, "--", "Estimated \mu", "LineWidth", 1.5);
    xlabel("\mu, with fitted variances held fixed");
    ylabel("Negative log-likelihood minus its fitted value");
    title("One-dimensional likelihood check near the optimum");
    grid on;
end

%% 12. Monte Carlo validation of mu estimator against CRLB

% The CRLB is a lower bound on the sampling VARIANCE of an unbiased
% estimator; it is not a bound on the squared error seen in one simulation.
% We therefore need many independent simulated experiments to estimate the
% repeated-sampling variance of the actual numerical ML estimator.
%
% Bias and variance answer different questions. Bias measures whether the
% estimator is centered on muTrue across experiments, whereas variance
% measures how widely its estimates fluctuate around their own mean.
%
% For this model, q = b'*b and
%
%   CRLB_mu = (tauTrue^2 + sigmaTrue^2/q)/T.
%
% Efficiency is CRLB_mu divided by the empirical variance. A value near one
% indicates that the estimator is operating close to the theoretical
% variance lower bound (subject to finite Monte Carlo sampling error).

rng(monteCarloSeed, "twister");

muHatMC = nan(nMonteCarlo, 1);
sigma2HatMC = nan(nMonteCarlo, 1);
exitFlagMC = nan(nMonteCarlo, 1);
successfulFitMC = false(nMonteCarlo, 1);

% Use the same optimizer settings as the single fit, except suppress output
% so the loop does not print one optimizer report per repetition.
optimizerOptionsMC = optimset(optimizerOptions, "Display", "off");

for r = 1:nMonteCarlo
    cTrueMC = muTrue + tauTrue * randn(T, 1);
    SMC = b * cTrueMC.' + sigmaTrue * randn(K, T);

    try
        % This is the same numerical ML fitting helper and likelihood used
        % for the original single-simulation analysis in Sections 5-6.
        fitMC = FitMarginalGaussianIID(SMC, b, optimizerOptionsMC);

        muHatMC(r) = fitMC.muHat;
        sigma2HatMC(r) = fitMC.sigma2Hat;
        exitFlagMC(r) = fitMC.exitFlag;

        successfulFitMC(r) = fitMC.exitFlag > 0 ...
            && all(isfinite(fitMC.thetaHat)) ...
            && all(isfinite([fitMC.muHat, fitMC.tau2Hat, fitMC.sigma2Hat])) ...
            && fitMC.tau2Hat > 0 && fitMC.sigma2Hat > 0 ...
            && isfinite(fitMC.negativeLogLikelihood);
    catch
        % A failed numerical fit is explicitly marked and excluded below.
        exitFlagMC(r) = -999;
    end
end

nSuccessfulFits = sum(successfulFitMC);
nFailedFits = nMonteCarlo - nSuccessfulFits;

if nSuccessfulFits < 2
    error("At least two successful Monte Carlo fits are required.");
end

muHatMCSuccessful = muHatMC(successfulFitMC);
meanMuHat = mean(muHatMCSuccessful);
biasMu = meanMuHat - muTrue;
varMuHat = var(muHatMCSuccessful, 0); % Sample normalization: 1/(R-1)
sdMuHat = std(muHatMCSuccessful, 0);

q = b.' * b;
crlbMu = (tauTrue^2 + sigmaTrue^2 / q) / T;
crlbSdMu = sqrt(crlbMu);
efficiencyMu = crlbMu / varMuHat;

zMu = (muHatMCSuccessful - muTrue) / crlbSdMu;
meanZMu = mean(zMu);
varZMu = var(zMu, 0);
mcSeMeanMu = sdMuHat / sqrt(nSuccessfulFits);

% Use exactly the same successful joint fits for sigma^2. The CRLB below is
% a bound on Var(sigma2Hat), not on Var(sqrt(sigma2Hat)) and not on the
% estimation error from any one simulated experiment.
sigma2HatMCSuccessful = sigma2HatMC(successfulFitMC);
trueSigma2 = sigmaTrue^2;
meanSigma2Hat = mean(sigma2HatMCSuccessful);
biasSigma2 = meanSigma2Hat - trueSigma2;
varSigma2Hat = var(sigma2HatMCSuccessful, 0); % Sample normalization: 1/(R-1)
sdSigma2Hat = std(sigma2HatMCSuccessful, 0);
mseSigma2 = mean((sigma2HatMCSuccessful - trueSigma2).^2);

% With tau^2 treated as an unknown nuisance parameter, sigma^2 is learned
% from the K-1 spectral directions orthogonal to b. Away from the tau^2 = 0
% boundary, T*(K-1)*sigma2Hat/sigmaTrue^2 has a chi-square distribution with
% T*(K-1) degrees of freedom. Thus sigma2Hat is unbiased and its variance is
% the following CRLB.
crlbSigma2 = 2 * sigmaTrue^4 / (T * (K - 1));
crlbSdSigma2 = sqrt(crlbSigma2);
efficiencySigma2 = crlbSigma2 / varSigma2Hat;

zSigma2 = (sigma2HatMCSuccessful - trueSigma2) / crlbSdSigma2;
meanZSigma2 = mean(zSigma2);
varZSigma2 = var(zSigma2, 0);
mcSeMeanSigma2 = sdSigma2Hat / sqrt(nSuccessfulFits);

fprintf("\nMonte Carlo validation of mu estimator against CRLB\n");
fprintf("Requested repetitions              = %d\n", nMonteCarlo);
fprintf("Successful fits                    = %d\n", nSuccessfulFits);
fprintf("Failed fits                        = %d\n", nFailedFits);
fprintf("True mu                            = %.8f\n", muTrue);
fprintf("Mean estimated mu                  = %.8f\n", meanMuHat);
fprintf("Estimated bias                     = %.8f\n", biasMu);
fprintf("Empirical variance of muHat        = %.8g\n", varMuHat);
fprintf("Empirical SD of muHat              = %.8f\n", sdMuHat);
fprintf("Analytical CRLB_mu                 = %.8g\n", crlbMu);
fprintf("sqrt(CRLB_mu)                      = %.8f\n", crlbSdMu);
fprintf("Efficiency (CRLB / empirical var) = %.8f\n", efficiencyMu);
fprintf("Mean normalized estimate zMu       = %.8f\n", meanZMu);
fprintf("Variance of zMu                    = %.8f\n", varZMu);
fprintf("Monte Carlo SE of mean(muHat)      = %.8f\n", mcSeMeanMu);

% Because sigma^2 is unbiased in this ideal model, comparing its empirical
% Monte Carlo variance with CRLB_sigma2 is a valid efficiency check. An
% efficiency near one means that the numerical joint ML fit extracts
% essentially all information available about sigma^2. A large absolute
% estimator variance alone would not imply poor performance if its CRLB
% were similarly large.
fprintf("\nMonte Carlo validation of sigma^2 estimator against CRLB\n");
fprintf("True sigma^2                       = %.8f\n", trueSigma2);
fprintf("Mean estimated sigma^2             = %.8f\n", meanSigma2Hat);
fprintf("Estimated bias                     = %.8f\n", biasSigma2);
fprintf("Empirical variance of sigma2Hat    = %.8g\n", varSigma2Hat);
fprintf("Empirical SD of sigma2Hat          = %.8f\n", sdSigma2Hat);
fprintf("Empirical MSE of sigma2Hat         = %.8g\n", mseSigma2);
fprintf("Analytical CRLB_sigma2             = %.8g\n", crlbSigma2);
fprintf("sqrt(CRLB_sigma2)                  = %.8f\n", crlbSdSigma2);
fprintf("Efficiency (CRLB / empirical var) = %.8f\n", efficiencySigma2);
fprintf("Mean normalized estimate zSigma2   = %.8f\n", meanZSigma2);
fprintf("Variance of zSigma2                = %.8f\n", varZSigma2);
fprintf("Monte Carlo SE of mean(sigma2Hat)  = %.8f\n", mcSeMeanSigma2);

% Histogram of the numerical ML estimates, with the Gaussian density
% predicted by N(muTrue, CRLB_mu). The density is written explicitly, so no
% Statistics and Machine Learning Toolbox function is required.
nHistogramBins = max(15, round(sqrt(nSuccessfulFits)));
muDensityGrid = linspace(muTrue - 4 * crlbSdMu, ...
    muTrue + 4 * crlbSdMu, 400);
predictedMuDensity = exp(-0.5 * ((muDensityGrid - muTrue) ...
    / crlbSdMu).^2) / (sqrt(2 * pi) * crlbSdMu);

figure("Name", "Monte Carlo distribution of muHat", "Color", "w");
histogram(muHatMCSuccessful, nHistogramBins, "Normalization", "pdf", ...
    "FaceColor", [0.35, 0.65, 0.85], "DisplayName", "ML estimates");
hold on;
plot(muDensityGrid, predictedMuDensity, "k-", "LineWidth", 2, ...
    "DisplayName", "N(\mu_{true}, CRLB_\mu)");
xline(muTrue, "r--", "\mu_{true}", "LineWidth", 1.7, ...
    "DisplayName", "True \mu");
xline(meanMuHat, "b-.", "Mean estimate", "LineWidth", 1.7, ...
    "DisplayName", "Mean estimated \mu");
xlabel("Numerical ML estimate \hat{\mu}");
ylabel("Probability density");
title("Monte Carlo distribution of \hat{\mu} (" ...
    + nSuccessfulFits + " successful fits)");
legend("Location", "best");
grid on;

% If the analytical result describes the numerical estimator well, these
% normalized estimates should resemble N(0,1).
zDensityGrid = linspace(-4, 4, 400);
standardNormalDensity = exp(-0.5 * zDensityGrid.^2) / sqrt(2 * pi);

figure("Name", "Normalized Monte Carlo estimates", "Color", "w");
histogram(zMu, nHistogramBins, "Normalization", "pdf", ...
    "FaceColor", [0.45, 0.75, 0.45], "DisplayName", "Normalized estimates");
hold on;
plot(zDensityGrid, standardNormalDensity, "k-", "LineWidth", 2, ...
    "DisplayName", "N(0,1)");
xline(0, "r--", "Zero", "LineWidth", 1.5, ...
    "DisplayName", "Expected mean");
xlabel("z_\mu = (\hat{\mu} - \mu_{true}) / sqrt(CRLB_\mu)");
ylabel("Probability density");
title("Normalized numerical ML estimates");
legend("Location", "best");
grid on;

% The finite-sample distribution of a variance estimator is not exactly
% Gaussian. Here T*(K-1) is large, so the Gaussian approximation based on
% the analytical CRLB should nevertheless be very accurate.
sigma2DensityGrid = linspace(trueSigma2 - 4 * crlbSdSigma2, ...
    trueSigma2 + 4 * crlbSdSigma2, 400);
predictedSigma2Density = exp(-0.5 * ((sigma2DensityGrid - trueSigma2) ...
    / crlbSdSigma2).^2) / (sqrt(2 * pi) * crlbSdSigma2);

figure("Name", "Monte Carlo distribution of sigma2Hat", "Color", "w");
histogram(sigma2HatMCSuccessful, nHistogramBins, ...
    "Normalization", "pdf", "FaceColor", [0.65, 0.50, 0.80], ...
    "DisplayName", "Joint ML estimates");
hold on;
plot(sigma2DensityGrid, predictedSigma2Density, "k-", ...
    "LineWidth", 2, "DisplayName", "Gaussian CRLB approximation");
xline(trueSigma2, "r--", "True \sigma^2", "LineWidth", 1.7, ...
    "DisplayName", "True \sigma^2");
xline(meanSigma2Hat, "b-.", "Mean estimate", "LineWidth", 1.7, ...
    "DisplayName", "Mean estimated \sigma^2");
xlabel("Numerical joint ML estimate \hat{\sigma}^2");
ylabel("Probability density");
title("Monte Carlo distribution of \hat{\sigma}^2 (" ...
    + nSuccessfulFits + " successful fits)");
legend("Location", "best");
grid on;

figure("Name", "Normalized Monte Carlo sigma2 estimates", "Color", "w");
histogram(zSigma2, nHistogramBins, "Normalization", "pdf", ...
    "FaceColor", [0.80, 0.60, 0.35], ...
    "DisplayName", "Normalized \sigma^2 estimates");
hold on;
plot(zDensityGrid, standardNormalDensity, "k-", "LineWidth", 2, ...
    "DisplayName", "N(0,1)");
xline(0, "r--", "Zero", "LineWidth", 1.5, ...
    "DisplayName", "Expected mean");
xlabel("z_{\sigma^2} = (\hat{\sigma}^2 - \sigma_{true}^2) / sqrt(CRLB_{\sigma^2})");
ylabel("Probability density");
title("Normalized numerical joint ML estimates of \sigma^2");
legend("Location", "best");
grid on;

%% Local fitting and likelihood functions

function fit = FitMarginalGaussianIID(S, b, optimizerOptions)
%FITMARGINALGAUSSIANIID Fit mu, tau^2, and sigma^2 with numerical ML.
%   This helper supplies identical initialization and likelihood logic to
%   the single experiment and every Monte Carlo repetition.

    likelihoodStatistics = BuildLikelihoodStatistics(S, b);
    q = likelihoodStatistics.q;
    cOLS = (likelihoodStatistics.bTS / q).';

    muInitial = mean(cOLS);
    orthogonalEnergy = likelihoodStatistics.sumSquaredS ...
        - sum(likelihoodStatistics.bTS.^2) / q;
    sigma2Initial = max(orthogonalEnergy ...
        / (likelihoodStatistics.K * likelihoodStatistics.T), 1e-6);
    tau2Initial = max(var(cOLS, 1) - sigma2Initial / q, 1e-6);
    thetaInitial = [muInitial, log(tau2Initial), log(sigma2Initial)];

    objective = @(theta) MarginalNegativeLogLikelihood( ...
        theta, likelihoodStatistics);
    [thetaHat, negativeLogLikelihood, exitFlag, optimizerOutput] = ...
        fminsearch(objective, thetaInitial, optimizerOptions);

    fit.thetaHat = thetaHat;
    fit.muHat = thetaHat(1);
    fit.tau2Hat = exp(thetaHat(2));
    fit.sigma2Hat = exp(thetaHat(3));
    fit.negativeLogLikelihood = negativeLogLikelihood;
    fit.exitFlag = exitFlag;
    fit.optimizerOutput = optimizerOutput;
    fit.thetaInitial = thetaInitial;
    fit.cOLS = cOLS;
    fit.likelihoodStatistics = likelihoodStatistics;
end

function statistics = BuildLikelihoodStatistics(S, b)
%BUILDLIKELIHOODSTATISTICS Precompute exact likelihood sufficient statistics.

    [statistics.K, statistics.T] = size(S);
    statistics.q = b.' * b;
    statistics.bTS = b.' * S;       % 1 x T projected measurements
    statistics.sumSquaredS = sum(S(:).^2);
end

function negativeLogLikelihood = MarginalNegativeLogLikelihood(theta, statistics)
%MARGINALNEGATIVELOGLIKELIHOOD Full IID marginal Gaussian objective.
%   The rank-one covariance V structure is used for speed. Along b, V has
%   eigenvalue sigma^2 + tau^2*(b'*b); in each of the K-1 orthogonal
%   directions it has eigenvalue sigma^2. This gives exactly the same log
%   determinant and quadratic form as a Cholesky solve of the full K x K V,
%   without forming inv(V) or repeatedly factoring V during Monte Carlo.

    K = statistics.K;
    T = statistics.T;
    q = statistics.q;
    bTS = statistics.bTS;

    mu = theta(1);
    tau2 = exp(theta(2));
    sigma2 = exp(theta(3));

    % Guard fminsearch against unusable trial points caused by exponentiation.
    if ~isfinite(mu) || ~isfinite(tau2) || ~isfinite(sigma2) || ...
            tau2 <= 0 || sigma2 <= 0
        negativeLogLikelihood = realmax("double") / 1e100;
        return;
    end

    parallelEigenvalue = sigma2 + tau2 * q;
    logDetV = (K - 1) * log(sigma2) + log(parallelEigenvalue);

    totalResidualEnergy = statistics.sumSquaredS ...
        - 2 * mu * sum(bTS) + T * mu^2 * q;
    parallelResidualEnergy = sum((bTS - mu * q).^2) / q;
    orthogonalResidualEnergy = max( ...
        totalResidualEnergy - parallelResidualEnergy, 0);

    quadraticTerm = orthogonalResidualEnergy / sigma2 ...
        + parallelResidualEnergy / parallelEigenvalue;

    negativeLogLikelihood = ...
        (T / 2) * logDetV ...
        + (1 / 2) * quadraticTerm ...
        + (T * K / 2) * log(2 * pi);
end
