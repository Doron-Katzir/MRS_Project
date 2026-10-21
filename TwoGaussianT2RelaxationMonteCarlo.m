%% TwoGaussianT2RelaxationMonteCarlo
% Standalone comparison of three T2 estimators for two overlapping Gaussians.
% Set T2_MC_TEST=1 before running for 8 repetitions at d=[0.2 0.4 3].
% T2_MC_FORCE_BASE=1 tests the fallback; T2_MC_OUTPUT_SUFFIX isolates output.
% Ordinary execution uses all 500 repetitions and 15 separations.
clear; clc; close all;

K = 160;
x = linspace(0, 10, K).';
width = 0.55;
dValues = linspace(0.2, 3, 15);
te = (10:20:190).';
sigma = 0.05;
nMC = 500;
seed = 1;
trueA = [1 1];
trueT2 = [100 100; 100 300];
scenarioNames = ["Equal T2", "Different T2"];
methodNames = ["Naive OLS NLS", "Corrected moment", "Joint NLS"];
if strcmp(getenv('T2_MC_TEST'), '1')
    nMC = 8;
    dValues = [0.2 0.4 3];
end
rng(seed, 'twister');
outDir = fullfile(fileparts(mfilename('fullpath')), ...
    ['T2RelaxationMonteCarloResults' getenv('T2_MC_OUTPUT_SUFFIX')]);
if ~isfolder(outDir), mkdir(outDir); end
useLsq = exist('lsqnonlin', 'file') == 2 && ...
    license('test', 'Optimization_Toolbox') && ...
    ~strcmp(getenv('T2_MC_FORCE_BASE'),'1');
if useLsq
    lsqOpts = optimoptions('lsqnonlin', 'Display', 'off', ...
        'FunctionTolerance', 1e-10, 'StepTolerance', 1e-10, ...
        'OptimalityTolerance', 1e-9, 'MaxIterations', 150, ...
        'MaxFunctionEvaluations', 1000);
else
    lsqOpts = [];
end

nD = numel(dValues); nSc = size(trueT2, 1); nTe = numel(te);
raw = struct('seed', seed, 'nMC', nMC, 'dValues', dValues, ...
    'te', te, 'sigma', sigma, 'trueA', trueA, 'trueT2', trueT2, ...
    'scenarioNames', scenarioNames, 'methodNames', methodNames, ...
    'solver', string(ternary(useLsq, 'lsqnonlin', 'bounded base MATLAB')));
% Dimension order: spectral point, echo, repetition, separation, scenario.
raw.spectra = zeros(K, nTe, nMC, nD, nSc);
raw.olsAmplitude = nan(2, nTe, nMC, nD, nSc);
% Dimension order: repetition, separation, scenario, metabolite, method.
raw.t2 = nan(nMC, nD, nSc, 2, 3);
raw.initialAmplitude = nan(nMC, nD, nSc, 2, 3);
raw.valid = false(nMC, nD, nSc, 2, 3);
raw.optimizerFailed = false(nMC, nD, nSc, 2, 3);
raw.boundHit = false(nMC, nD, nSc, 2, 3);
raw.invalidMoment = false(nMC, nD, nSc, 2);
raw.momentQ = nan(nMC, nD, nSc, 2);
raw.momentDenominator = nan(nMC, nD, nSc, 2);
raw.noiseCovariance = nan(2, 2, nD);

maxOlsError = 0; maxCovarianceError = 0;
for id = 1:nD
    B = [exp(-0.5*((x-(5-dValues(id)/2))/width).^2), ...
         exp(-0.5*((x-(5+dValues(id)/2))/width).^2)];
    % B\S uses QR rather than forming the ill-conditioned normal equations.
    L = B \ eye(K);
    V = sigma^2 * (L * L.');
    analyticV = sigma^2 * ((B.'*B) \ eye(2));
    raw.noiseCovariance(:, :, id) = V;
    maxCovarianceError = max(maxCovarianceError, ...
        norm(V-analyticV, 'fro'));
    for isc = 1:nSc
        aTrue = exp(-te ./ trueT2(isc, :)) .* trueA;
        cleanS = B * aTrue.';
        maxOlsError = max(maxOlsError, ...
            max(abs(B\cleanS-aTrue.'), [], 'all'));
        for imc = 1:nMC
            S = cleanS + sigma*randn(K, nTe);
            raw.spectra(:, :, imc, id, isc) = S;
            aHat = B \ S;
            raw.olsAmplitude(:, :, imc, id, isc) = aHat;
            separateFit = nan(2, 2);
            for j = 1:2
                [p, ok, hit] = fitExponential(te, aHat(j, :).', ...
                    useLsq, lsqOpts);
                separateFit(j, :) = p;
                raw.initialAmplitude(imc, id, isc, j, 1) = p(1);
                raw.t2(imc, id, isc, j, 1) = p(2);
                raw.optimizerFailed(imc, id, isc, j, 1) = ~ok;
                raw.boundHit(imc, id, isc, j, 1) = hit;
                raw.valid(imc, id, isc, j, 1) = ok && ~hit;

                v = V(j, j);
                y = aHat(j, :).';
                denom = sum(y(1:end-1).^2) - (nTe-1)*v;
                q = sum(y(1:end-1).*y(2:end)) / denom;
                raw.momentDenominator(imc, id, isc, j) = denom;
                raw.momentQ(imc, id, isc, j) = q;
                validMoment = isfinite(q) && denom > 0 && q > 0 && q < 1;
                raw.invalidMoment(imc, id, isc, j) = ~validMoment;
                if validMoment
                    raw.t2(imc, id, isc, j, 2) = ...
                        -(te(2)-te(1))/log(q);
                    raw.valid(imc, id, isc, j, 2) = true;
                end
            end
            [pJoint, ok, hit] = fitJoint(te, S, B, separateFit, ...
                useLsq, lsqOpts);
            for j = 1:2
                raw.initialAmplitude(imc, id, isc, j, 3) = pJoint(j);
                raw.t2(imc, id, isc, j, 3) = pJoint(j+2);
                raw.optimizerFailed(imc, id, isc, j, 3) = ~ok;
                raw.boundHit(imc, id, isc, j, 3) = hit(j);
                raw.valid(imc, id, isc, j, 3) = ok && ~hit(j);
            end
        end
        fprintf('Experiment 1: scenario %d, d=%.1f, %d repetitions complete\n', ...
            isc, dValues(id), nMC);
    end
end

assert(maxOlsError < 1e-10, 'Noise-free OLS validation failed.');
assert(maxCovarianceError < 1e-10, 'OLS noise covariance validation failed.');
validation = table(maxOlsError, maxCovarianceError, useLsq, nMC, seed, ...
    'VariableNames', {'MaxNoiseFreeOLSError', 'MaxNoiseCovarianceError', ...
    'UsedLsqnonlin', 'Repetitions', 'Seed'});
writetable(validation, fullfile(outDir, 'Validation.csv'));

rowTemplate = struct('Scenario', "", 'Separation', NaN, 'Metabolite', 0, ...
    'Method', "", 'Bias_ms', NaN, 'RMSE_ms', NaN, ...
    'MedianAbsError_ms', NaN, 'ValidPercent', NaN, 'ValidCount', 0, ...
    'OptimizerFailureCount', 0, 'BoundHitCount', 0, ...
    'InvalidMomentCount', 0);
rows = repmat(rowTemplate, nSc*nD*2*3, 1);
rowIndex = 0;
metric = nan(nSc, nD, 3, 2, 4); % bias, RMSE, median AE, valid percent
for isc = 1:nSc
    for id = 1:nD
        for j = 1:2
            for m = 1:3
                v = squeeze(raw.valid(:, id, isc, j, m));
                estimates = squeeze(raw.t2(:, id, isc, j, m));
                errors = estimates(v)-trueT2(isc, j);
                bias = NaN; rmse = NaN; medae = NaN;
                if ~isempty(errors)
                    bias = mean(errors);
                    rmse = sqrt(mean(errors.^2));
                    medae = median(abs(errors));
                end
                validPct = 100*sum(v)/nMC;
                metric(isc, id, m, j, :) = [bias rmse medae validPct];
                rowIndex = rowIndex+1;
                rows(rowIndex) = struct('Scenario', scenarioNames(isc), ...
                    'Separation', dValues(id), 'Metabolite', j, ...
                    'Method', methodNames(m), 'Bias_ms', bias, ...
                    'RMSE_ms', rmse, 'MedianAbsError_ms', medae, ...
                    'ValidPercent', validPct, 'ValidCount', sum(v), ...
                    'OptimizerFailureCount', sum(raw.optimizerFailed(:,id,isc,j,m)), ...
                    'BoundHitCount', sum(raw.boundHit(:,id,isc,j,m)), ...
                    'InvalidMomentCount', sum(raw.invalidMoment(:,id,isc,j))*(m==2));
            end
        end
    end
end
summary = struct2table(rows);
writetable(summary, fullfile(outDir, 'T2Metrics.csv'));
raw.validation = validation;
raw.metric = metric;
save(fullfile(outDir, 'RawResults.mat'), 'raw', 'summary', '-v7.3');

plotMetric(metric, dValues, scenarioNames, methodNames, 2, ...
    'T2 RMSE (ms)', 'T2_RMSE', outDir);
plotMetric(metric, dValues, scenarioNames, methodNames, 4, ...
    'Valid estimates (%)', 'ValidRate', outDir);
fprintf('Solver: %s; repetitions: %d; seed: %d\n', raw.solver, nMC, seed);
fprintf('Maximum noise-free OLS error: %.3e\n', maxOlsError);
fprintf('Noise covariance identity error: %.3e\n', maxCovarianceError);
fprintf('Saved experiment 1 to %s\n', outDir);

function [p, ok, hit] = fitExponential(te, y, useLsq, opts)
    lb = [1e-8 10]; ub = [20 2000];
    p0 = [max(0.1, min(5, max(y(1), mean(y)))), 180];
    if useLsq
        [p, ~, ~, exitflag] = lsqnonlin( ...
            @(z) z(1)*exp(-te/z(2))-y, p0, lb, ub, opts);
    else
        % Profile out amplitude, leaving a bounded one-dimensional NLS fit.
        f = @(t) profileSSE(t, te, y, lb(1), ub(1));
        [t, ~, exitflag] = fminbnd(f, lb(2), ub(2), ...
            optimset('Display','off','TolX',1e-9,'MaxFunEvals',1000));
        g = exp(-te/t);
        p = [min(ub(1),max(lb(1),(g.'*y)/(g.'*g))), t];
    end
    ok = exitflag > 0 && all(isfinite(p));
    hit = any(p <= lb+1e-5.*(ub-lb) | p >= ub-1e-5.*(ub-lb));
end

function sse = profileSSE(t, te, y, lo, hi)
    g = exp(-te/t);
    a = min(hi, max(lo, (g.'*y)/(g.'*g)));
    r = a*g-y;
    sse = r.'*r;
end

function [p, ok, hit] = fitJoint(te, S, B, separateFit, useLsq, opts)
    lb = [1e-8 1e-8 10 10]; ub = [20 20 2000 2000];
    p0 = [0.8 0.8 180 180];
    if all(isfinite(separateFit), 'all')
        p0 = [separateFit(:,1).' separateFit(:,2).'];
        p0 = min(ub-1e-4, max(lb+1e-4, p0));
    end
    objective = @(z) reshape(B*(exp(-te./z(3:4)).*z(1:2)).'-S, [], 1);
    if useLsq
        [p, ~, ~, exitflag] = lsqnonlin(objective, p0, lb, ub, opts);
    else
        p0 = min(ub-1e-7, max(lb+1e-7, p0));
        u0 = log((p0-lb)./(ub-p0));
        wrapped = @(u) sum(objective(lb+(ub-lb)./(1+exp(-u))).^2);
        [u, ~, exitflag] = fminsearch(wrapped, u0, ...
            optimset('Display','off','TolX',1e-9,'TolFun',1e-10, ...
            'MaxFunEvals',3000,'MaxIter',1000));
        p = lb+(ub-lb)./(1+exp(-u));
    end
    ok = exitflag > 0 && all(isfinite(p));
    hit = [any(p([1 3]) <= lb([1 3])+1e-5.*(ub([1 3])-lb([1 3])) | ...
               p([1 3]) >= ub([1 3])-1e-5.*(ub([1 3])-lb([1 3]))), ...
           any(p([2 4]) <= lb([2 4])+1e-5.*(ub([2 4])-lb([2 4])) | ...
               p([2 4]) >= ub([2 4])-1e-5.*(ub([2 4])-lb([2 4])) )];
end

function plotMetric(metric, dValues, scenarios, methods, index, yText, fileStem, outDir)
    colors = [0 0.447 0.741; 0.850 0.325 0.098; 0.470 0.200 0.620];
    fig = figure('Visible','off','Color','w','Position',[100 100 1100 720]);
    tiledlayout(2,2,'TileSpacing','compact','Padding','compact');
    for isc = 1:2
        for j = 1:2
            ax = nexttile;
            set(ax,'Color','w','XColor','k','YColor','k');
            hold on;
            for m = 1:3
                plot(dValues, squeeze(metric(isc,:,m,j,index)), '-o', ...
                    'Color',colors(m,:),'LineWidth',1.4,'MarkerSize',4, ...
                    'DisplayName',methods(m));
            end
            grid on;
            xlabel('Gaussian separation d','Color','k');
            ylabel(yText,'Color','k');
            title(sprintf('%s, metabolite %d',scenarios(isc),j),'Color','k');
            if isc==1 && j==1
                legend('Location','best','Color','w','TextColor','k');
            end
            if index==4, ylim([0 105]); end
        end
    end
    exportgraphics(fig, fullfile(outDir, [fileStem '.png']), 'Resolution', 180);
    savefig(fig, fullfile(outDir, [fileStem '.fig']));
    close(fig);
end

function out = ternary(condition, a, b)
    if condition, out = a; else, out = b; end
end
