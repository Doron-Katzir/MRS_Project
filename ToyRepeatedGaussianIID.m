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
tau2HatMC = nan(nMonteCarlo, 1);
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
        tau2HatMC(r) = fitMC.tau2Hat;
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

% Use the same successful joint fits to study tau^2. Unlike mu and sigma^2,
% the interior ML estimator of tau^2 is downward biased at finite T. We
% therefore compare its empirical bias, variance, and MSE with their
% finite-sample predictions rather than report an ordinary CRLB efficiency.
tau2HatMCSuccessful = tau2HatMC(successfulFitMC);
trueTau2 = tauTrue^2;
meanTau2Hat = mean(tau2HatMCSuccessful);
biasTau2 = meanTau2Hat - trueTau2;
varTau2Hat = var(tau2HatMCSuccessful, 0); % Sample normalization: 1/(R-1)
sdTau2Hat = std(tau2HatMCSuccessful, 0);
mseTau2 = mean((tau2HatMCSuccessful - trueTau2).^2);
mcSeMeanTau2 = sdTau2Hat / sqrt(nSuccessfulFits);

vTrue = tauTrue^2 + sigmaTrue^2 / q;
biasTau2Analytic = -vTrue / T;
expectedTau2HatAnalytic = trueTau2 + biasTau2Analytic;
varTau2Analytic = 2 * (T - 1) / T^2 * vTrue^2 ...
    + 2 * sigmaTrue^4 / (T * (K - 1) * q^2);
mseTau2Analytic = varTau2Analytic + biasTau2Analytic^2;
meanErrorVsPrediction = meanTau2Hat - expectedTau2HatAnalytic;

biasRatio = biasTau2 / biasTau2Analytic;
varianceRatio = varTau2Hat / varTau2Analytic;
mseRatio = mseTau2 / mseTau2Analytic;

% The log parameterization keeps tau^2 strictly positive. Values at or below
% this documented numerical threshold are classified as boundary-near. The
% interior finite-sample formulas cease to be exact if such estimates occur
% frequently.
tau2BoundaryThreshold = 1e-6;
isTau2BoundaryNear = tau2HatMCSuccessful <= tau2BoundaryThreshold;
nTau2BoundaryNear = sum(isTau2BoundaryNear);
percentTau2BoundaryNear = 100 * nTau2BoundaryNear / nSuccessfulFits;

% This is the ordinary CRLB for UNBIASED estimators of tau^2. Since the ML
% estimator studied here is biased at finite T, CRLB_tau2/varTau2Hat is not
% an appropriate efficiency measure. The bound is printed only as a
% theoretical variance reference while bias, variance, and MSE are assessed.
crlbTau2Unbiased = (2 / T) * (vTrue^2 ...
    + sigmaTrue^4 / (q^2 * (K - 1)));

% The predicted bias magnitude is vTrue/T, so it decreases as 1/T. A later
% experiment could check this at T = 36, 50, 100, and 200; no additional
% simulations or fits are performed here.

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

% Agreement of the empirical and analytical bias, variance, and MSE shows
% whether the numerical joint ML estimator follows finite-sample theory.
% Ratios near one indicate agreement; biasRatio should not be overinterpreted
% in settings where the predicted bias is extremely close to zero.
fprintf("\nFinite-sample validation of tau^2 ML estimator\n");
fprintf("True tau^2                              = %.8f\n", trueTau2);
fprintf("Monte Carlo mean tau^2 estimate         = %.8f\n", meanTau2Hat);
fprintf("Analytical expected tau^2 estimate      = %.8f\n", ...
    expectedTau2HatAnalytic);
fprintf("Mean error versus analytical prediction = %.8f\n", ...
    meanErrorVsPrediction);
fprintf("Empirical bias                          = %.8f\n", biasTau2);
fprintf("Analytical bias                         = %.8f\n", ...
    biasTau2Analytic);
fprintf("Bias ratio (empirical / analytical)     = %.8f\n", biasRatio);
fprintf("Empirical variance                      = %.8g\n", varTau2Hat);
fprintf("Analytical variance                     = %.8g\n", ...
    varTau2Analytic);
fprintf("Variance ratio (empirical / analytical) = %.8f\n", ...
    varianceRatio);
fprintf("Empirical MSE                           = %.8g\n", mseTau2);
fprintf("Analytical MSE                          = %.8g\n", ...
    mseTau2Analytic);
fprintf("MSE ratio (empirical / analytical)      = %.8f\n", mseRatio);
fprintf("Empirical SD                            = %.8f\n", sdTau2Hat);
fprintf("Monte Carlo SE of mean(tau2Hat)         = %.8f\n", ...
    mcSeMeanTau2);
fprintf("Ordinary CRLB for unbiased estimators   = %.8g\n", ...
    crlbTau2Unbiased);
fprintf("Successful fits                         = %d\n", nSuccessfulFits);
fprintf("Failed fits                             = %d\n", nFailedFits);
fprintf("Boundary-near threshold                 = %.1e\n", ...
    tau2BoundaryThreshold);
fprintf("Boundary-near tau^2 estimates           = %d (%.4f%%)\n", ...
    nTau2BoundaryNear, percentTau2BoundaryNear);

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

% The tau^2 histogram displays the finite-sample downward shift directly.
figure("Name", "Monte Carlo distribution of tau2Hat", "Color", "w");
histogram(tau2HatMCSuccessful, nHistogramBins, ...
    "FaceColor", [0.35, 0.70, 0.70], "DisplayName", "Joint ML estimates");
hold on;
xline(trueTau2, "r--", "True \tau^2", "LineWidth", 1.7, ...
    "DisplayName", "True \tau^2");
xline(meanTau2Hat, "b-.", "Monte Carlo mean", "LineWidth", 1.7, ...
    "DisplayName", "Monte Carlo mean");
xline(expectedTau2HatAnalytic, "k:", "Analytical E[\hat{\tau}^2]", ...
    "LineWidth", 2, "DisplayName", "Analytical expected value");
xlabel("Numerical joint ML estimate \hat{\tau}^2");
ylabel("Number of Monte Carlo estimates");
title("Finite-sample distribution of \hat{\tau}^2");
legend("Location", "best");
grid on;

% The error histogram compares the observed mean error with the predicted
% finite-sample bias. Their agreement is more relevant here than comparison
% of the biased estimator variance with the ordinary unbiased CRLB.
tau2EstimationError = tau2HatMCSuccessful - trueTau2;

figure("Name", "Monte Carlo tau2 estimation error", "Color", "w");
histogram(tau2EstimationError, nHistogramBins, ...
    "FaceColor", [0.85, 0.55, 0.35], "DisplayName", "Estimation errors");
hold on;
xline(0, "r--", "Zero error", "LineWidth", 1.7, ...
    "DisplayName", "Zero error");
xline(biasTau2, "b-.", "Empirical mean error", "LineWidth", 1.7, ...
    "DisplayName", "Empirical bias");
xline(biasTau2Analytic, "k:", "Analytical bias", "LineWidth", 2, ...
    "DisplayName", "Analytical bias");
xlabel("\hat{\tau}^2 - \tau_{true}^2");
ylabel("Number of Monte Carlo estimates");
title("Finite-sample estimation error of \hat{\tau}^2");
legend("Location", "best");
grid on;

%% 13. Effect of number of repeated measurements T

% This experiment combines analytical finite-sample predictions with Monte
% Carlo fits at several values of T. Increasing T supplies more independent
% realizations of c_t and therefore more information about the population
% parameters. The basis, K, and all true parameters remain unchanged.
TValues = [12, 18, 24, 36, 50, 72, 100, 200];
nMonteCarloT = 500;
monteCarloTSeed = 456;

TGrid = TValues(:);
nTValues = numel(TGrid);
qT = b.' * b;
vTrueT = tauTrue^2 + sigmaTrue^2 / qT;
trueTau2T = tauTrue^2;
trueSigma2T = sigmaTrue^2;

% Analytical predictions. The variances of mu and sigma^2 decrease as 1/T.
% Noise estimation is especially precise because each spectrum supplies K-1
% directions orthogonal to b. In contrast, tau^2 describes between-time
% concentration variability and is fundamentally learned from only T draws.
muBiasAnalyticT = zeros(nTValues, 1);
muVarianceAnalyticT = vTrueT ./ TGrid;
muMSEAnalyticT = muVarianceAnalyticT;

sigma2BiasAnalyticT = zeros(nTValues, 1);
sigma2VarianceAnalyticT = 2 * sigmaTrue^4 ./ (TGrid * (K - 1));
sigma2MSEAnalyticT = sigma2VarianceAnalyticT;

% The interior ML bias of tau^2 falls as 1/T, while its standard deviation
% falls approximately as 1/sqrt(T).
tau2BiasAnalyticT = -vTrueT ./ TGrid;
tau2ExpectedAnalyticT = trueTau2T + tau2BiasAnalyticT;
tau2VarianceAnalyticT = 2 * (TGrid - 1) ./ TGrid.^2 * vTrueT^2 ...
    + 2 * sigmaTrue^4 ./ (TGrid * (K - 1) * qT^2);
tau2MSEAnalyticT = tau2VarianceAnalyticT + tau2BiasAnalyticT.^2;
crlbTau2UnbiasedT = (2 ./ TGrid) .* (vTrueT^2 ...
    + sigmaTrue^4 / (qT^2 * (K - 1)));

tau2RelativeBiasAnalyticT = abs(tau2BiasAnalyticT) / trueTau2T;
tau2RelativeSDAnalyticT = sqrt(tau2VarianceAnalyticT) / trueTau2T;
tau2RelativeRMSEAnalyticT = sqrt(tau2MSEAnalyticT) / trueTau2T;

% Preallocate empirical summaries, one entry for each value of T.
muBiasEmpiricalT = nan(nTValues, 1);
muVarianceEmpiricalT = nan(nTValues, 1);
muMSEEmpiricalT = nan(nTValues, 1);

sigma2BiasEmpiricalT = nan(nTValues, 1);
sigma2VarianceEmpiricalT = nan(nTValues, 1);
sigma2MSEEmpiricalT = nan(nTValues, 1);

tau2BiasEmpiricalT = nan(nTValues, 1);
tau2VarianceEmpiricalT = nan(nTValues, 1);
tau2MSEEmpiricalT = nan(nTValues, 1);
tau2RelativeBiasEmpiricalT = nan(nTValues, 1);
tau2RelativeSDEmpiricalT = nan(nTValues, 1);
tau2RelativeRMSEEmpiricalT = nan(nTValues, 1);

nSuccessfulFitsT = zeros(nTValues, 1);
nFailedFitsT = zeros(nTValues, 1);
nTau2BoundaryNearT = zeros(nTValues, 1);

rng(monteCarloTSeed, "twister");
optimizerOptionsT = optimset(optimizerOptions, "Display", "off");

for tIndex = 1:nTValues
    TCurrent = TGrid(tIndex);

    muEstimateCurrent = nan(nMonteCarloT, 1);
    tau2EstimateCurrent = nan(nMonteCarloT, 1);
    sigma2EstimateCurrent = nan(nMonteCarloT, 1);
    successfulCurrent = false(nMonteCarloT, 1);

    for r = 1:nMonteCarloT
        cTrueCurrent = muTrue + tauTrue * randn(TCurrent, 1);
        SCurrent = b * cTrueCurrent.' ...
            + sigmaTrue * randn(K, TCurrent);

        try
            % Use the same numerical joint ML estimator as Sections 5-6 and
            % 12; no analytical shortcut replaces the fitted estimates.
            fitCurrent = FitMarginalGaussianIID( ...
                SCurrent, b, optimizerOptionsT);

            muEstimateCurrent(r) = fitCurrent.muHat;
            tau2EstimateCurrent(r) = fitCurrent.tau2Hat;
            sigma2EstimateCurrent(r) = fitCurrent.sigma2Hat;

            successfulCurrent(r) = fitCurrent.exitFlag > 0 ...
                && all(isfinite(fitCurrent.thetaHat)) ...
                && all(isfinite([fitCurrent.muHat, ...
                    fitCurrent.tau2Hat, fitCurrent.sigma2Hat])) ...
                && fitCurrent.tau2Hat > 0 && fitCurrent.sigma2Hat > 0 ...
                && isfinite(fitCurrent.negativeLogLikelihood);
        catch
            successfulCurrent(r) = false;
        end
    end

    nSuccessfulFitsT(tIndex) = sum(successfulCurrent);
    nFailedFitsT(tIndex) = nMonteCarloT - nSuccessfulFitsT(tIndex);

    if nSuccessfulFitsT(tIndex) < 2
        warning("Fewer than two successful fits for T = %d.", TCurrent);
        continue;
    end

    muSuccessful = muEstimateCurrent(successfulCurrent);
    tau2Successful = tau2EstimateCurrent(successfulCurrent);
    sigma2Successful = sigma2EstimateCurrent(successfulCurrent);

    muBiasEmpiricalT(tIndex) = mean(muSuccessful) - muTrue;
    muVarianceEmpiricalT(tIndex) = var(muSuccessful, 0);
    muMSEEmpiricalT(tIndex) = mean((muSuccessful - muTrue).^2);

    sigma2BiasEmpiricalT(tIndex) = mean(sigma2Successful) - trueSigma2T;
    sigma2VarianceEmpiricalT(tIndex) = var(sigma2Successful, 0);
    sigma2MSEEmpiricalT(tIndex) = mean( ...
        (sigma2Successful - trueSigma2T).^2);

    tau2BiasEmpiricalT(tIndex) = mean(tau2Successful) - trueTau2T;
    tau2VarianceEmpiricalT(tIndex) = var(tau2Successful, 0);
    tau2MSEEmpiricalT(tIndex) = mean( ...
        (tau2Successful - trueTau2T).^2);
    tau2RelativeBiasEmpiricalT(tIndex) = ...
        abs(tau2BiasEmpiricalT(tIndex)) / trueTau2T;
    tau2RelativeSDEmpiricalT(tIndex) = ...
        sqrt(tau2VarianceEmpiricalT(tIndex)) / trueTau2T;
    tau2RelativeRMSEEmpiricalT(tIndex) = ...
        sqrt(tau2MSEEmpiricalT(tIndex)) / trueTau2T;

    nTau2BoundaryNearT(tIndex) = sum( ...
        tau2Successful <= tau2BoundaryThreshold);
end

% One row per T summarizes both the Monte Carlo results and their analytical
% targets. CRLB_tau2_unbiased is a reference only, not an efficiency measure
% for the biased finite-sample ML estimator.
tScalingSummary = table( ...
    TGrid, ...
    muBiasEmpiricalT, muVarianceEmpiricalT, muVarianceAnalyticT, ...
    muMSEEmpiricalT, ...
    sigma2BiasEmpiricalT, sigma2VarianceEmpiricalT, ...
    sigma2VarianceAnalyticT, sigma2MSEEmpiricalT, ...
    tau2BiasEmpiricalT, tau2BiasAnalyticT, ...
    tau2VarianceEmpiricalT, tau2VarianceAnalyticT, ...
    tau2MSEEmpiricalT, tau2MSEAnalyticT, crlbTau2UnbiasedT, ...
    tau2RelativeBiasEmpiricalT, tau2RelativeBiasAnalyticT, ...
    tau2RelativeSDEmpiricalT, tau2RelativeSDAnalyticT, ...
    tau2RelativeRMSEEmpiricalT, tau2RelativeRMSEAnalyticT, ...
    nTau2BoundaryNearT, nSuccessfulFitsT, nFailedFitsT, ...
    VariableNames=[ ...
        "T", ...
        "muBias_MC", "muVariance_MC", "muVariance_CRLB", "muMSE_MC", ...
        "sigma2Bias_MC", "sigma2Variance_MC", "sigma2Variance_CRLB", ...
        "sigma2MSE_MC", ...
        "tau2Bias_MC", "tau2Bias_Analytic", ...
        "tau2Variance_MC", "tau2Variance_Analytic", ...
        "tau2MSE_MC", "tau2MSE_Analytic", "tau2CRLB_Unbiased", ...
        "tau2RelativeBias_MC", "tau2RelativeBias_Analytic", ...
        "tau2RelativeSD_MC", "tau2RelativeSD_Analytic", ...
        "tau2RelativeRMSE_MC", "tau2RelativeRMSE_Analytic", ...
        "nTau2BoundaryNear", "nSuccessfulFits", "nFailedFits"]);

fprintf("\nEffect of number of repeated measurements T\n");
disp(tScalingSummary);

% Figure 1: mu variance and its CRLB.
figure("Name", "Effect of T on mu variance", "Color", "w");
plot(TGrid, muVarianceAnalyticT, "k-o", "LineWidth", 1.8, ...
    "DisplayName", "Analytical CRLB");
hold on;
plot(TGrid, muVarianceEmpiricalT, "b-s", "LineWidth", 1.5, ...
    "DisplayName", "Monte Carlo variance");
xlabel("Number of repeated spectra T");
ylabel("Variance of \hat{\mu}");
title("Precision of population-mean estimation versus T");
legend("Location", "best");
grid on;

% Figure 2: sigma^2 variance and its CRLB.
figure("Name", "Effect of T on sigma2 variance", "Color", "w");
plot(TGrid, sigma2VarianceAnalyticT, "k-o", "LineWidth", 1.8, ...
    "DisplayName", "Analytical CRLB");
hold on;
plot(TGrid, sigma2VarianceEmpiricalT, "b-s", "LineWidth", 1.5, ...
    "DisplayName", "Monte Carlo variance");
xlabel("Number of repeated spectra T");
ylabel("Variance of \hat{\sigma}^2");
title("Precision of spectral-noise variance estimation versus T");
legend("Location", "best");
grid on;

% Figure 3: the approximately 1/T finite-sample bias of tau^2.
figure("Name", "Effect of T on tau2 bias", "Color", "w");
plot(TGrid, tau2BiasAnalyticT, "k-o", "LineWidth", 1.8, ...
    "DisplayName", "Analytical bias");
hold on;
plot(TGrid, tau2BiasEmpiricalT, "b-s", "LineWidth", 1.5, ...
    "DisplayName", "Monte Carlo bias");
yline(0, "r--", "Zero bias", "LineWidth", 1.2, ...
    "DisplayName", "Zero");
xlabel("Number of repeated spectra T");
ylabel("Bias of \hat{\tau}^2");
title("Finite-sample ML bias of \tau^2 versus T");
legend("Location", "best");
grid on;

% Figure 4: tau^2 variance. The ordinary unbiased CRLB is shown only as a
% reference and is not used to define efficiency for the biased ML estimator.
figure("Name", "Effect of T on tau2 variance", "Color", "w");
plot(TGrid, tau2VarianceAnalyticT, "k-o", "LineWidth", 1.8, ...
    "DisplayName", "Analytical ML variance");
hold on;
plot(TGrid, tau2VarianceEmpiricalT, "b-s", "LineWidth", 1.5, ...
    "DisplayName", "Monte Carlo ML variance");
plot(TGrid, crlbTau2UnbiasedT, "r--^", "LineWidth", 1.4, ...
    "DisplayName", "Ordinary unbiased-estimator CRLB");
xlabel("Number of repeated spectra T");
ylabel("Variance of \hat{\tau}^2");
title("Population-variance estimation uncertainty versus T");
legend("Location", "best");
grid on;

% Figure 5: tau^2 MSE combines finite-sample bias and variance.
figure("Name", "Effect of T on tau2 MSE", "Color", "w");
plot(TGrid, tau2MSEAnalyticT, "k-o", "LineWidth", 1.8, ...
    "DisplayName", "Analytical MSE");
hold on;
plot(TGrid, tau2MSEEmpiricalT, "b-s", "LineWidth", 1.5, ...
    "DisplayName", "Monte Carlo MSE");
xlabel("Number of repeated spectra T");
ylabel("MSE of \hat{\tau}^2");
title("Mean-squared error of \tau^2 estimation versus T");
legend("Location", "best");
grid on;

% Figure 6: relative RMSE directly quantifies uncertainty compared with the
% true population variance. T=36 is highlighted because it is the repeated
% scan count in the current MRS dataset; adequacy depends on required precision.
figure("Name", "Relative uncertainty of tau2 versus T", "Color", "w");
plot(TGrid, 100 * tau2RelativeRMSEAnalyticT, "k-o", ...
    "LineWidth", 1.8, "DisplayName", "Analytical relative RMSE");
hold on;
plot(TGrid, 100 * tau2RelativeRMSEEmpiricalT, "b-s", ...
    "LineWidth", 1.5, "DisplayName", "Monte Carlo relative RMSE");
xline(36, "r--", "T = 36", "LineWidth", 1.5, ...
    "DisplayName", "Current MRS T");
xlabel("Number of repeated spectra T");
ylabel("RMSE relative to true \tau^2 (%)");
title("Relative uncertainty of population-variance estimation");
legend("Location", "best");
grid on;

% Focused quantitative report for the current MRS repeat count. No binary
% claim is made about whether 36 is enough; these values allow adequacy to be
% judged against the precision required by the scientific application.
t36Index = find(TGrid == 36, 1);
if ~isempty(t36Index)
    fprintf("\nFocused results for T = 36\n");
    fprintf("MU:\n");
    fprintf("  Monte Carlo SD             = %.8f\n", ...
        sqrt(muVarianceEmpiricalT(t36Index)));
    fprintf("  Analytical CRLB SD         = %.8f\n", ...
        sqrt(muVarianceAnalyticT(t36Index)));
    fprintf("SIGMA^2:\n");
    fprintf("  Monte Carlo SD             = %.8f\n", ...
        sqrt(sigma2VarianceEmpiricalT(t36Index)));
    fprintf("  Analytical CRLB SD         = %.8f\n", ...
        sqrt(sigma2VarianceAnalyticT(t36Index)));
    fprintf("TAU^2:\n");
    fprintf("  Monte Carlo bias           = %.8f\n", ...
        tau2BiasEmpiricalT(t36Index));
    fprintf("  Analytical expected bias   = %.8f\n", ...
        tau2BiasAnalyticT(t36Index));
    fprintf("  Monte Carlo SD             = %.8f\n", ...
        sqrt(tau2VarianceEmpiricalT(t36Index)));
    fprintf("  Analytical SD              = %.8f\n", ...
        sqrt(tau2VarianceAnalyticT(t36Index)));
    fprintf("  Monte Carlo RMSE           = %.8f\n", ...
        sqrt(tau2MSEEmpiricalT(t36Index)));
    fprintf("  Analytical RMSE            = %.8f\n", ...
        sqrt(tau2MSEAnalyticT(t36Index)));
    fprintf("  Monte Carlo relative bias  = %.2f%% of true tau^2\n", ...
        100 * tau2RelativeBiasEmpiricalT(t36Index));
    fprintf("  Analytical relative bias   = %.2f%% of true tau^2\n", ...
        100 * tau2RelativeBiasAnalyticT(t36Index));
    fprintf("  Monte Carlo relative SD    = %.2f%% of true tau^2\n", ...
        100 * tau2RelativeSDEmpiricalT(t36Index));
    fprintf("  Analytical relative SD     = %.2f%% of true tau^2\n", ...
        100 * tau2RelativeSDAnalyticT(t36Index));
    fprintf("  Monte Carlo relative RMSE  = %.2f%% of true tau^2\n", ...
        100 * tau2RelativeRMSEEmpiricalT(t36Index));
    fprintf("  Analytical relative RMSE   = %.2f%% of true tau^2\n", ...
        100 * tau2RelativeRMSEAnalyticT(t36Index));
end

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
