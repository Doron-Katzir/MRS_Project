%% TwoGaussianT2CovarianceMonteCarlo
% Standalone covariance-subtraction experiment for two Gaussian T2 decays.
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
t2MethodNames = ["Uncorrected covariance", "Corrected covariance", ...
    "Conventional NLS"];
covMethodNames = ["Uncorrected OLS", "Corrected OLS", "Model-based NLS"];
if strcmp(getenv('T2_MC_TEST'), '1')
    nMC = 8;
    dValues = [0.2 0.4 3];
end
rng(seed, 'twister');
outDir = fullfile(fileparts(mfilename('fullpath')), ...
    ['T2CovarianceMonteCarloResults' getenv('T2_MC_OUTPUT_SUFFIX')]);
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
    'scenarioNames', scenarioNames, 't2MethodNames', t2MethodNames, ...
    'covMethodNames', covMethodNames, ...
    'solver', string(ternary(useLsq, 'lsqnonlin', 'bounded base MATLAB')));
% Dimension order: spectral point, echo, repetition, separation, scenario.
raw.spectra = zeros(K, nTe, nMC, nD, nSc);
raw.olsAmplitude = nan(2, nTe, nMC, nD, nSc);
raw.fittedAmplitude = nan(2, nTe, nMC, nD, nSc);
raw.meanAmplitude = nan(2, nMC, nD, nSc);
raw.noiseCovariance = nan(2, 2, nD);
raw.trueCovariance = nan(2, 2, nSc);
% Covariance: row, column, repetition, separation, scenario, method.
raw.covariance = nan(2, 2, nMC, nD, nSc, 3);
raw.frobeniusError = nan(nMC, nD, nSc, 3);
raw.negativeEigenvalue = false(nMC, nD, nSc);
raw.negativeDiagonal = false(nMC, nD, nSc, 2);
% T2: repetition, separation, scenario, metabolite, method.
raw.t2 = nan(nMC, nD, nSc, 2, 3);
raw.initialAmplitude = nan(nMC, nD, nSc, 2);
raw.ratio = nan(nMC, nD, nSc, 2, 2);
raw.valid = false(nMC, nD, nSc, 2, 3);
raw.invalidRatio = false(nMC, nD, nSc, 2, 2);
raw.optimizerFailed = false(nMC, nD, nSc, 2);
raw.boundHit = false(nMC, nD, nSc, 2);

ratioGrid = logspace(log10(10), log10(2000), 2001);
ratioValues = ratioFunction(ratioGrid, te);
ratioMonotone = all(diff(ratioValues) < 0);
assert(ratioMonotone, 'R(T) must decrease on [10, 2000] ms.');
maxRatioRecoveryError = 0;
for isc = 1:nSc
    aTrue = exp(-te ./ trueT2(isc, :)) .* trueA;
    raw.trueCovariance(:, :, isc) = cov(aTrue);
    for j = 1:2
        [recovered, valid] = invertRatio( ...
            ratioFunction(trueT2(isc,j), te), te);
        assert(valid, 'True T2 is outside the ratio inversion interval.');
        maxRatioRecoveryError = max(maxRatioRecoveryError, ...
            abs(recovered-trueT2(isc,j)));
    end
end
assert(maxRatioRecoveryError < 1e-6, 'Ratio inversion validation failed.');

maxOlsError = 0; maxCovarianceError = 0;
for id = 1:nD
    B = [exp(-0.5*((x-(5-dValues(id)/2))/width).^2), ...
         exp(-0.5*((x-(5+dValues(id)/2))/width).^2)];
    L = B \ eye(K);
    V = sigma^2 * (L*L.');
    analyticV = sigma^2 * ((B.'*B) \ eye(2));
    raw.noiseCovariance(:, :, id) = V;
    maxCovarianceError = max(maxCovarianceError, ...
        norm(V-analyticV, 'fro'));
    for isc = 1:nSc
        aTrue = exp(-te ./ trueT2(isc, :)) .* trueA;
        cleanS = B*aTrue.';
        maxOlsError = max(maxOlsError, ...
            max(abs(B\cleanS-aTrue.'), [], 'all'));
        sTrue = raw.trueCovariance(:, :, isc);
        for imc = 1:nMC
            S = cleanS + sigma*randn(K, nTe);
            raw.spectra(:, :, imc, id, isc) = S;
            aHat = B\S;
            raw.olsAmplitude(:, :, imc, id, isc) = aHat;
            mu = mean(aHat, 2);
            raw.meanAmplitude(:, imc, id, isc) = mu;
            sOLS = cov(aHat.');
            sCorrected = sOLS-V;
            raw.negativeEigenvalue(imc, id, isc) = any(eig(sCorrected)<0);
            raw.negativeDiagonal(imc, id, isc, :) = diag(sCorrected)<0;
            aFit = nan(nTe, 2);
            for j = 1:2
                [p, ok, hit] = fitExponential(te, aHat(j,:).', ...
                    useLsq, lsqOpts);
                raw.initialAmplitude(imc, id, isc, j) = p(1);
                raw.t2(imc, id, isc, j, 3) = p(2);
                raw.optimizerFailed(imc, id, isc, j) = ~ok;
                raw.boundHit(imc, id, isc, j) = hit;
                raw.valid(imc, id, isc, j, 3) = ok && ~hit;
                if ok
                    aFit(:,j) = p(1)*exp(-te/p(2));
                end
                for m = 1:2
                    if m==1, s = sOLS; else, s = sCorrected; end
                    rHat = s(j,j)/(mu(j)^2);
                    raw.ratio(imc, id, isc, j, m) = rHat;
                    [tHat, isValid] = invertRatio(rHat, te);
                    raw.t2(imc, id, isc, j, m) = tHat;
                    raw.valid(imc, id, isc, j, m) = isValid;
                    raw.invalidRatio(imc, id, isc, j, m) = ~isValid;
                end
            end
            raw.fittedAmplitude(:, :, imc, id, isc) = aFit.';
            if all(isfinite(aFit), 'all')
                sNLS = cov(aFit);
            else
                sNLS = nan(2);
            end
            sAll = cat(3, sOLS, sCorrected, sNLS);
            for m = 1:3
                raw.covariance(:, :, imc, id, isc, m) = sAll(:,:,m);
                raw.frobeniusError(imc,id,isc,m) = ...
                    norm(sAll(:,:,m)-sTrue, 'fro');
            end
        end
        fprintf('Experiment 2: scenario %d, d=%.1f, %d repetitions complete\n', ...
            isc, dValues(id), nMC);
    end
end
assert(maxOlsError < 1e-10, 'Noise-free OLS validation failed.');
assert(maxCovarianceError < 1e-10, 'OLS noise covariance validation failed.');

% Compare Monte Carlo covariance means with their exact expectations.
% Five empirical standard errors allow for simultaneous checks across cells.
expectationTemplate = struct('Scenario',"",'Separation',NaN,'Method',"", ...
    'MatrixRow',0,'MatrixColumn',0,'ObservedMean',NaN, ...
    'ExpectedMean',NaN,'StandardError',NaN,'ZScore',NaN,'Pass',false);
expectationRows = repmat(expectationTemplate,nSc*nD*2*4,1);
expectationIndex = 0;
meanCovariance = nan(2,2,nD,nSc,3);
covarianceBias = nan(2,2,nD,nSc,3);
for isc = 1:nSc
    for id = 1:nD
        for m = 1:3
            c = raw.covariance(:,:,:,id,isc,m);
            c = reshape(c,2,2,nMC);
            meanCovariance(:,:,id,isc,m) = mean(c,3,'omitnan');
            covarianceBias(:,:,id,isc,m) = ...
                meanCovariance(:,:,id,isc,m)-raw.trueCovariance(:,:,isc);
            if m>2, continue; end
            target = raw.trueCovariance(:,:,isc);
            if m==1, target = target+raw.noiseCovariance(:,:,id); end
            for row = 1:2
                for col = 1:2
                    samples = squeeze(c(row,col,:));
                    se = std(samples)/sqrt(nMC);
                    observed = mean(samples);
                    z = abs(observed-target(row,col))/max(se,eps);
                    expectationIndex = expectationIndex+1;
                    expectationRows(expectationIndex) = struct( ...
                        'Scenario',scenarioNames(isc), ...
                        'Separation',dValues(id), ...
                        'Method',covMethodNames(m), ...
                        'MatrixRow',row,'MatrixColumn',col, ...
                        'ObservedMean',observed, ...
                        'ExpectedMean',target(row,col), ...
                        'StandardError',se,'ZScore',z, ...
                        'Pass',z<=5);
                end
            end
        end
    end
end
expectationChecks = struct2table(expectationRows);
writetable(expectationChecks, fullfile(outDir,'CovarianceExpectationChecks.csv'));
raw.meanCovariance = meanCovariance;
raw.covarianceBias = covarianceBias;

t2Template = struct('Scenario',"",'Separation',NaN,'Metabolite',0, ...
    'Method',"",'Bias_ms',NaN,'RMSE_ms',NaN,'MedianAbsError_ms',NaN, ...
    'ValidPercent',NaN,'ValidCount',0,'InvalidRatioCount',0, ...
    'OptimizerFailureCount',0,'BoundHitCount',0);
t2Rows = repmat(t2Template,nSc*nD*2*3,1);
t2Index = 0;
t2Metric = nan(nSc,nD,3,2,4); % bias, RMSE, median AE, valid percent
for isc = 1:nSc
    for id = 1:nD
        for j = 1:2
            for m = 1:3
                valid = squeeze(raw.valid(:,id,isc,j,m));
                estimates = squeeze(raw.t2(:,id,isc,j,m));
                errors = estimates(valid)-trueT2(isc,j);
                bias = NaN; rmse = NaN; medae = NaN;
                if ~isempty(errors)
                    bias = mean(errors);
                    rmse = sqrt(mean(errors.^2));
                    medae = median(abs(errors));
                end
                validPct = 100*sum(valid)/nMC;
                t2Metric(isc,id,m,j,:) = [bias rmse medae validPct];
                invalidCount = 0; failureCount = 0; hitCount = 0;
                if m<=2
                    invalidCount = sum(raw.invalidRatio(:,id,isc,j,m));
                else
                    failureCount = sum(raw.optimizerFailed(:,id,isc,j));
                    hitCount = sum(raw.boundHit(:,id,isc,j));
                end
                t2Index = t2Index+1;
                t2Rows(t2Index) = struct('Scenario',scenarioNames(isc), ...
                    'Separation',dValues(id),'Metabolite',j, ...
                    'Method',t2MethodNames(m),'Bias_ms',bias, ...
                    'RMSE_ms',rmse,'MedianAbsError_ms',medae, ...
                    'ValidPercent',validPct,'ValidCount',sum(valid), ...
                    'InvalidRatioCount',invalidCount, ...
                    'OptimizerFailureCount',failureCount, ...
                    'BoundHitCount',hitCount);
            end
        end
    end
end
t2Summary = struct2table(t2Rows);
writetable(t2Summary, fullfile(outDir,'T2Metrics.csv'));

covTemplate = struct('Scenario',"",'Separation',NaN,'Method',"", ...
    'RMSE_F',NaN,'ValidCount',0,'NegativeEigenvalueCount',0, ...
    'NegativeDiagonal1Count',0,'NegativeDiagonal2Count',0);
covRows = repmat(covTemplate,nSc*nD*3,1);
covIndex = 0;
covMetric = nan(nSc,nD,3);
matrixTemplate = struct('Scenario',"",'Separation',NaN,'Method',"", ...
    'MatrixRow',0,'MatrixColumn',0,'TrueValue',NaN, ...
    'MeanEstimate',NaN,'Bias',NaN);
matrixRows = repmat(matrixTemplate,nSc*nD*3*4,1);
matrixIndex = 0;
for isc = 1:nSc
    for id = 1:nD
        for m = 1:3
            e = raw.frobeniusError(:,id,isc,m);
            finite = isfinite(e);
            covMetric(isc,id,m) = sqrt(mean(e(finite).^2));
            negativeEigenCount = 0; negativeDiag = [0 0];
            if m==2
                negativeEigenCount = sum(raw.negativeEigenvalue(:,id,isc));
                negativeDiag = [sum(raw.negativeDiagonal(:,id,isc,1)), ...
                    sum(raw.negativeDiagonal(:,id,isc,2))];
            end
            covIndex = covIndex+1;
            covRows(covIndex) = struct('Scenario',scenarioNames(isc), ...
                'Separation',dValues(id),'Method',covMethodNames(m), ...
                'RMSE_F',covMetric(isc,id,m),'ValidCount',sum(finite), ...
                'NegativeEigenvalueCount',negativeEigenCount, ...
                'NegativeDiagonal1Count',negativeDiag(1), ...
                'NegativeDiagonal2Count',negativeDiag(2));
            for row = 1:2
                for col = 1:2
                    matrixIndex = matrixIndex+1;
                    matrixRows(matrixIndex) = struct( ...
                        'Scenario',scenarioNames(isc), ...
                        'Separation',dValues(id), ...
                        'Method',covMethodNames(m), ...
                        'MatrixRow',row,'MatrixColumn',col, ...
                        'TrueValue',raw.trueCovariance(row,col,isc), ...
                        'MeanEstimate',meanCovariance(row,col,id,isc,m), ...
                        'Bias',covarianceBias(row,col,id,isc,m));
                end
            end
        end
    end
end
covSummary = struct2table(covRows);
matrixSummary = struct2table(matrixRows);
writetable(covSummary, fullfile(outDir,'CovarianceMetrics.csv'));
writetable(matrixSummary, fullfile(outDir,'CovarianceMatrixBias.csv'));

validation = table(maxOlsError,maxCovarianceError,maxRatioRecoveryError, ...
    ratioMonotone,sum(expectationChecks.Pass),height(expectationChecks), ...
    max(expectationChecks.ZScore),useLsq,nMC,seed, ...
    'VariableNames',{'MaxNoiseFreeOLSError','MaxNoiseCovarianceError', ...
    'MaxRatioRecoveryError_ms','RatioMonotone', ...
    'ExpectationPassCount','ExpectationCheckCount','MaxExpectationZ', ...
    'UsedLsqnonlin','Repetitions','Seed'});
writetable(validation, fullfile(outDir,'Validation.csv'));
raw.validation = validation;
raw.t2Metric = t2Metric;
raw.covMetric = covMetric;
save(fullfile(outDir,'RawResults.mat'),'raw','t2Summary', ...
    'covSummary','matrixSummary','expectationChecks','-v7.3');

plotT2Metric(t2Metric,dValues,scenarioNames,t2MethodNames,2, ...
    'T2 RMSE (ms)','T2_RMSE',outDir);
plotT2Metric(t2Metric,dValues,scenarioNames,t2MethodNames,1, ...
    'T2 bias (ms)','T2_Bias',outDir);
plotCovMetric(covMetric,dValues,scenarioNames,covMethodNames,outDir);
plotMatrixComparison(raw,dValues,meanCovariance,scenarioNames,outDir);
fprintf('Solver: %s; repetitions: %d; seed: %d\n',raw.solver,nMC,seed);
fprintf('Maximum noise-free OLS error: %.3e\n',maxOlsError);
fprintf('Noise covariance identity error: %.3e\n',maxCovarianceError);
fprintf('Ratio recovery error: %.3e ms; monotone: %d\n', ...
    maxRatioRecoveryError,ratioMonotone);
fprintf('Covariance expectation checks passed: %d/%d (max |z| %.2f)\n', ...
    sum(expectationChecks.Pass),height(expectationChecks), ...
    max(expectationChecks.ZScore));
fprintf('Corrected covariance negative eigenvalue counts: %d, %d\n', ...
    sum(raw.negativeEigenvalue(:,:,1),'all'), ...
    sum(raw.negativeEigenvalue(:,:,2),'all'));
fprintf('Saved experiment 2 to %s\n',outDir);

function r = ratioFunction(t, te)
    g = exp(-te(:)./t(:).');
    mu = mean(g,1);
    v = sum((g-mu).^2,1)/(numel(te)-1);
    r = reshape(v./mu.^2,size(t));
end

function [tHat, valid] = invertRatio(rHat, te)
    tHat = NaN;
    rHigh = ratioFunction(10,te);
    rLow = ratioFunction(2000,te);
    valid = isfinite(rHat) && rHat>rLow && rHat<rHigh;
    if valid
        try
            tHat = fzero(@(t) ratioFunction(t,te)-rHat,[10 2000]);
            valid = isfinite(tHat) && tHat>=10 && tHat<=2000;
        catch
            valid = false;
        end
    end
    if ~valid, tHat = NaN; end
end

function [p, ok, hit] = fitExponential(te, y, useLsq, opts)
    lb = [1e-8 10]; ub = [20 2000];
    p0 = [max(0.1,min(5,max(y(1),mean(y)))),180];
    if useLsq
        [p,~,~,exitflag] = lsqnonlin( ...
            @(z) z(1)*exp(-te/z(2))-y,p0,lb,ub,opts);
    else
        f = @(t) profileSSE(t,te,y,lb(1),ub(1));
        [t,~,exitflag] = fminbnd(f,lb(2),ub(2), ...
            optimset('Display','off','TolX',1e-9,'MaxFunEvals',1000));
        g = exp(-te/t);
        p = [min(ub(1),max(lb(1),(g.'*y)/(g.'*g))),t];
    end
    ok = exitflag>0 && all(isfinite(p));
    hit = any(p<=lb+1e-5.*(ub-lb) | p>=ub-1e-5.*(ub-lb));
end

function sse = profileSSE(t,te,y,lo,hi)
    g = exp(-te/t);
    a = min(hi,max(lo,(g.'*y)/(g.'*g)));
    r = a*g-y;
    sse = r.'*r;
end

function plotT2Metric(metric,dValues,scenarios,methods,index,yText,stem,outDir)
    colors = [0 0.447 0.741; 0.850 0.325 0.098; 0.470 0.200 0.620];
    fig = figure('Visible','off','Color','w','Position',[100 100 1100 720]);
    tiledlayout(2,2,'TileSpacing','compact','Padding','compact');
    for isc = 1:2
        for j = 1:2
            ax = nexttile; hold on;
            set(ax,'Color','w','XColor','k','YColor','k');
            if index==1, yline(0,'k:','HandleVisibility','off'); end
            for m = 1:3
                plot(dValues,squeeze(metric(isc,:,m,j,index)),'-o', ...
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
        end
    end
    exportgraphics(fig,fullfile(outDir,[stem '.png']),'Resolution',180);
    savefig(fig,fullfile(outDir,[stem '.fig'])); close(fig);
end

function plotCovMetric(metric,dValues,scenarios,methods,outDir)
    colors = [0 0.447 0.741; 0.850 0.325 0.098; 0.470 0.200 0.620];
    fig = figure('Visible','off','Color','w','Position',[100 100 1100 430]);
    tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
    for isc = 1:2
        ax = nexttile; hold on;
        set(ax,'Color','w','XColor','k','YColor','k');
        for m = 1:3
            plot(dValues,squeeze(metric(isc,:,m)),'-o', ...
                'Color',colors(m,:),'LineWidth',1.4,'MarkerSize',4, ...
                'DisplayName',methods(m));
        end
        grid on; xlabel('Gaussian separation d','Color','k');
        ylabel('Covariance Frobenius RMSE','Color','k');
        title(scenarios(isc),'Color','k');
        if isc==1
            legend('Location','best','Color','w','TextColor','k');
        end
    end
    exportgraphics(fig,fullfile(outDir,'Covariance_RMSE.png'),'Resolution',180);
    savefig(fig,fullfile(outDir,'Covariance_RMSE.fig')); close(fig);
end

function plotMatrixComparison(raw,dValues,meanCov,scenarios,outDir)
    [~,id] = min(abs(dValues-0.4));
    fig = figure('Visible','off','Color','w','Position',[100 100 1400 650]);
    tiledlayout(2,4,'TileSpacing','compact','Padding','compact');
    labels = ["True covariance","Mean uncorrected", ...
        "Mean corrected","Mean model-based NLS"];
    for isc = 1:2
        matrices = cat(3,raw.trueCovariance(:,:,isc), ...
            meanCov(:,:,id,isc,1),meanCov(:,:,id,isc,2), ...
            meanCov(:,:,id,isc,3));
        lo = min(matrices,[],'all'); hi = max(matrices,[],'all');
        if hi==lo, hi=lo+eps; end
        for m = 1:4
            ax = nexttile; imagesc(matrices(:,:,m));
            set(ax,'Color','w','XColor','k','YColor','k');
            clim([lo hi]); axis image;
            xticks(1:2); yticks(1:2);
            xlabel('Metabolite','Color','k');
            ylabel('Metabolite','Color','k');
            title(sprintf('%s: %s',scenarios(isc),labels(m)),'Color','k');
            for row = 1:2
                for col = 1:2
                    value = matrices(row,col,m);
                    if (value-lo)/(hi-lo)>0.55
                        color = 'k';
                    else
                        color = 'w';
                    end
                    text(col,row,sprintf('%.5f',value), ...
                        'HorizontalAlignment','center','Color',color, ...
                        'FontWeight','bold','FontSize',9);
                end
            end
            cb = colorbar;
            cb.Color = 'k';
        end
    end
    sgtitle(sprintf('Covariance comparison at d = %.1f',dValues(id)), ...
        'Color','k');
    exportgraphics(fig,fullfile(outDir,'Covariance_Matrices.png'), ...
        'Resolution',180);
    savefig(fig,fullfile(outDir,'Covariance_Matrices.fig')); close(fig);
end

function out = ternary(condition,a,b)
    if condition, out=a; else, out=b; end
end
