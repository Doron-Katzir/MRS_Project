%% 
clear;
clc;

cfg = ProjectConfig();

%% User input for this run

cfg.paths.rootDir = "C:\Users\User\Documents\Thesis_Lab";
cfg.paths.dataDir = fullfile(cfg.paths.rootDir, "Data");
cfg.paths.coordDir = fullfile(cfg.paths.rootDir, "LCMFit");

%% Splice data and save coord files
%% Splice only one file

% cfg.input.mode = "singleFile";
% cfg.input.singlePatientID = "P01";
% cfg.input.singleFileName = "meas_MID00020_FID54986_eja_svs_slaser_TE_80_1_PRE.dat";
% cfg.input.singleFile = fullfile(cfg.paths.dataDir, cfg.input.singleFileName);

%% Splice all files in the working directory

% cfg.input.mode = "directory";
% 
% cfg.input.directory = cfg.paths.dataDir;
% cfg.input.filePattern = "*.dat";
% cfg.input.recursive = false;

%% Call splice function
% processedPatients = Splice_data_multi_patient_input(cfg);

%% Load coord files
%% Load only one fitted patient folder

% cfg.load.mode = "singleSubfolder";
% cfg.load.selectedPatientID = "P11";
% cfg.load.coordDir = cfg.paths.coordDir;
% 

%% Load all subfolders raw, coord & basis files

cfg.load.mode = "allSubfolders";
cfg.load.coordDir = cfg.paths.coordDir;
cfg.load.includeRawBasis = true;
cfg.load.loadOnly = true;

outputs = Load_subsets_multi_patient_input(cfg);


%% Examine signals
% Select patients and extract their data
patientIDs = ["P01"];
patientClasses = getPatients(outputs, patientIDs);

% Get FIDs, spectrums and basis matrices
[fidSigs, spectrumSigs, sigLengths, basisMat] = extractSignals( ...
    patientClasses);

%% Plotting example
patientID = "P01";  % Choose patient
metab = "NAA";  % Choose metab

plotFIDSpec(outputs, patientID, metab)


%% Construct ppm axis from acquisition metadata

[sFit, BFit, ppmFit] = getAxis( ...
    patientClasses, spectrumSigs, sigLengths, basisMat);

%% Complex nonnegative least squares
[cHat, sHat, residual] = nnlsFit(BFit, sFit, ppmFit, Plot=true);


%% Shared Lorentzian broadening experiment

dt = patient.preFitBasis.globalHeader.badelt;
n = size(B, 1);
t = (0:n-1).' * dt;

% Recover time-domain basis functions
basisFids = ifft(ifftshift(B, 1), [], 1);

% Candidate linewidths in Hz
gammaGrid = 0:0.5:30;

errors = zeros(size(gammaGrid));
fits = cell(size(gammaGrid));
concentrations = cell(size(gammaGrid));

y = [real(sFit); imag(sFit)];

for k = 1:numel(gammaGrid)

    gamma = gammaGrid(k);

    % Apply Lorentzian broadening in time
    decay = exp(-pi * gamma * t);
    broadenedFids = basisFids .* decay;

    % Convert to frequency domain
    Bgamma = fftshift(fft(broadenedFids, [], 1), 1);
    BgFit = Bgamma(mask, :);

    % Normalize basis columns
    scale = sqrt(sum(abs(BgFit).^2, 1));
    Bnorm = BgFit ./ scale;

    % Solve complex NNLS as a real problem
    A = [real(Bnorm); imag(Bnorm)];
    aHat = lsqnonneg(A, y);

    cHat = aHat ./ scale(:);
    sHat = BgFit * cHat;

    % Store results
    errors(k) = norm(sFit - sHat) / norm(sFit);
    fits{k} = sHat;
    concentrations{k} = cHat;

end

% Find optimal candidate
[bestError, bestIdx] = min(errors);

bestGamma = gammaGrid(bestIdx);
bestFit = fits{bestIdx};
bestC = concentrations{bestIdx};

fprintf('Original error: %.2f%%\n', 100*errors(1));
fprintf('Best gamma: %.2f Hz\n', bestGamma);
fprintf('Best error: %.2f%%\n', 100*bestError);


figure;
plot(gammaGrid, 100*errors, 'o-', 'LineWidth', 1.5);

xlabel('Additional Lorentzian linewidth (Hz)');
ylabel('Relative residual (%)');
title('Profiled NNLS Residual vs Lorentzian Broadening');
grid on;


%% Compare original and broadened fits

figure('Name', 'Effect of Lorentzian Broadening');
tiledlayout(3,1);

% Real component
nexttile;
plot(ppmFit, real(sFit), 'b', 'LineWidth', 1.2); hold on;
plot(ppmFit, real(fits{1}), 'r--');
plot(ppmFit, real(bestFit), 'g', 'LineWidth', 1.2);
set(gca, 'XDir', 'reverse');
legend('Measured', 'Original NNLS', 'Broadened NNLS');
title(sprintf('Real Spectrum | gamma = %.1f Hz', bestGamma));
grid on;

% Imaginary component
nexttile;
plot(ppmFit, imag(sFit), 'b', 'LineWidth', 1.2); hold on;
plot(ppmFit, imag(fits{1}), 'r--');
plot(ppmFit, imag(bestFit), 'g', 'LineWidth', 1.2);
set(gca, 'XDir', 'reverse');
legend('Measured', 'Original NNLS', 'Broadened NNLS');
title('Imaginary Spectrum');
grid on;

% Remaining residual
nexttile;
r = sFit - bestFit;
plot(ppmFit, real(r), 'b'); hold on;
plot(ppmFit, imag(r), 'r');
set(gca, 'XDir', 'reverse');
legend('Real residual', 'Imaginary residual');
xlabel('Chemical shift (ppm)');
title(sprintf('Remaining Residual: %.1f%%', 100*bestError));
grid on;



%% Global phase investigation

gamma = bestGamma;

% Construct the broadened basis
decay = exp(-pi * gamma * t);
Bbroadened = fftshift(fft(basisFids .* decay, [], 1), 1);
Bg = Bbroadened(mask, :);

% Normalize columns
scale = sqrt(sum(abs(Bg).^2, 1));
BgNorm = Bg ./ scale;

y = [real(sFit); imag(sFit)];

% Candidate global phases (degrees)
phaseGrid = -180:2:180;
phaseErrors = zeros(size(phaseGrid));
phaseFits = cell(size(phaseGrid));

for k = 1:numel(phaseGrid)

    phi = deg2rad(phaseGrid(k));

    % Rotate the complex basis
    Brot = exp(1i * phi) * BgNorm;

    % Re-estimate concentrations
    A = [real(Brot); imag(Brot)];
    aHat = lsqnonneg(A, y);

    cHat = aHat ./ scale(:);

    % Reconstruct
    sHat = exp(1i * phi) * Bg * cHat;

    phaseFits{k} = sHat;
    phaseErrors(k) = norm(sFit - sHat) / norm(sFit);

end

[bestPhaseError, idxPhase] = min(phaseErrors);
bestPhase = phaseGrid(idxPhase);
bestPhaseFit = phaseFits{idxPhase};

fprintf("Best phase: %.1f degrees\n", bestPhase);
fprintf("Residual without phase correction: %.2f%%\n", ...
    100 * phaseErrors(phaseGrid == 0));
fprintf("Residual with phase correction: %.2f%%\n", ...
    100 * bestPhaseError);

figure;
plot(phaseGrid, 100 * phaseErrors, 'LineWidth', 1.5);
xlabel("Global phase (degrees)");
ylabel("Relative residual (%)");
title("Residual vs Global Phase");
grid on;


%% Joint linewidth and phase grid search

gammaGrid2 = 0:2:24;       % Hz
phaseGrid2 = -80:5:20;     % degrees

errorMap = zeros(numel(gammaGrid2), numel(phaseGrid2));

y = [real(sFit); imag(sFit)];

for ig = 1:numel(gammaGrid2)

    gamma = gammaGrid2(ig);

    % Apply Lorentzian broadening
    decay = exp(-pi * gamma * t);
    Bgamma = fftshift(fft(basisFids .* decay, [], 1), 1);
    Bg = Bgamma(mask, :);

    % Normalize columns
    scale = sqrt(sum(abs(Bg).^2, 1));
    Bnorm = Bg ./ scale;

    for ip = 1:numel(phaseGrid2)

        phi = deg2rad(phaseGrid2(ip));

        % Global phase rotation
        Brot = exp(1i * phi) * Bnorm;

        % NNLS
        A = [real(Brot); imag(Brot)];
        aHat = lsqnonneg(A, y);

        % Relative residual
        errorMap(ig, ip) = norm(y - A*aHat) / norm(y);

    end
end

% Joint optimum on the grid
[bestErrorJoint, linearIdx] = min(errorMap(:));
[ig, ip] = ind2sub(size(errorMap), linearIdx);

gammaJoint = gammaGrid2(ig);
phaseJoint = phaseGrid2(ip);

fprintf('Best linewidth: %.1f Hz\n', gammaJoint);
fprintf('Best phase: %.1f degrees\n', phaseJoint);
fprintf('Joint residual: %.2f%%\n', 100*bestErrorJoint);

% Visualize the objective surface

figure;
contourf(phaseGrid2, gammaGrid2, 100*errorMap, 20);
colorbar;
hold on;
plot(phaseJoint, gammaJoint, 'rx', ...
    'MarkerSize', 12, 'LineWidth', 2);

xlabel('Global phase (degrees)');
ylabel('Lorentzian linewidth (Hz)');
title('Joint Profiled NNLS Residual (%)');


%% Reconstruct the jointly optimized model

decay = exp(-pi * gammaJoint * t);

Bgamma = fftshift(fft(basisFids .* decay, [], 1), 1);
Bg = Bgamma(mask, :);

phi = deg2rad(phaseJoint);
Bjoint = exp(1i * phi) * Bg;

scale = vecnorm(Bjoint, 2, 1);
Bnorm = Bjoint ./ scale;

A = [real(Bnorm); imag(Bnorm)];
y = [real(sFit); imag(sFit)];

aHat = lsqnonneg(A, y);
cJoint = aHat ./ scale(:);

sJoint = Bjoint * cJoint;
rJoint = sFit - sJoint;

% Plot

figure;
tiledlayout(3,1);

nexttile;
plot(ppmFit, real(sFit), 'b'); hold on;
plot(ppmFit, real(sJoint), 'r');
set(gca, 'XDir', 'reverse');
legend('Measured', 'Joint fit');
title('Real Spectrum');
grid on;

nexttile;
plot(ppmFit, imag(sFit), 'b'); hold on;
plot(ppmFit, imag(sJoint), 'r');
set(gca, 'XDir', 'reverse');
legend('Measured', 'Joint fit');
title('Imaginary Spectrum');
grid on;

nexttile;
plot(ppmFit, real(rJoint), 'b'); hold on;
plot(ppmFit, imag(rJoint), 'r');
set(gca, 'XDir', 'reverse');
legend('Real residual', 'Imaginary residual');
title('Remaining Residual');
xlabel('Chemical shift (ppm)');
grid on;


%% Full all shifts, all broadening, global shift joint fit

%% ============================================================
%  NONLINEAR MRS SPECTRAL FITTING
%  Variable projection: lsqnonlin + lsqnonneg
%  ============================================================

%% 1. Select measurement from the cached dataset

patientIDs = string({outputs.patientResults.patientID});
p = find(patientIDs == "P01", 1);
assert(~isempty(p), "Patient not found.");

patient = outputs.patientResults(p);
T = patient.originalData.fileTable;

idx = find(T.division == 36 & T.part == 1, 1);
assert(~isempty(idx), "Measurement not found.");

fid = patient.originalData.fids(:, idx);
n = numel(fid);

s = fftshift(fft(fid));

%% 2. Recover original time-domain basis functions

basisStored = patient.preFitBasis.spectra;
metaboliteNames = string(patient.preFitBasis.metaboliteNames);
M = numel(metaboliteNames);

basisFids = complex(zeros(n, M));
tailRatios = zeros(M, 1);

for j = 1:M

    if iscell(basisStored)
        bj = basisStored{j}(:);
    else
        bj = basisStored(:, j);
    end

    % Undo fftshift and FFT of the exported basis
    bPadded = ifft(ifftshift(bj));

    assert(numel(bPadded) == 2*n, ...
        "Unexpected basis length.");

    % Verify that the second half is zero padding
    tailRatios(j) = norm(bPadded(n+1:end)) / ...
                    max(norm(bPadded(1:n)), eps);

    basisFids(:, j) = bPadded(1:n);
end

assert(max(tailRatios) < 1e-5, ...
    "Basis does not appear to be correctly zero-filled.");

%% 3. Construct frequency axis and fitting window

hdr = patient.preFitBasis.globalHeader;

dt = hdr.badelt;
t = (0:n-1).' * dt;

SW = 1/dt;
freqHz = -SW/2 + (0:n-1).' * SW/n;

B0 = hdr.hzpppm / SpinsJ.GetGyromagneticRatio("1H");

ppm = SpinsJ.HzToPPM( ...
    -freqHz, B0, "1H", SpinsJ.waterPPM);

mask = ppm >= 0.2 & ppm <= 4.0;

ppmFit = ppm(mask);
sFit = s(mask);

assert(norm(sFit) > 0, "Measured spectrum is empty.");

%% 4. Nonlinear parameter initialization

% theta = [gamma(1:M); delta(1:M); g; phi]
%
% gamma: Lorentzian FWHM in Hz
% delta: metabolite shifts in Hz
% g:     Gaussian FWHM in Hz
% phi:   global phase in radians

gammaInitial = 10;
gamma0 = gammaInitial * ones(M, 1);
gammaUpper = 80 * ones(M, 1);
gammaLower = 0 * ones(M, 1);
gammaScale = abs(gammaUpper - gammaLower);

shiftInitial = 0;
delta0 = shiftInitial * ones(M, 1);
shiftUpper = 25 * ones(M, 1);
shiftLower = -15 * ones(M, 1);
shiftScale = abs(shiftUpper - shiftLower);

globalBroad0 = 0;
globalBroadUpper = 25;
globalBroadLower = 0;
globalBroadScale = abs(globalBroadUpper - globalBroadLower);

phi0 = deg2rad(-35);
phiUpper = pi;
phiLower = -pi;
phiScale = abs(phiUpper - phiLower);

theta = [gamma0; delta0; globalBroad0; phi0];

% Bounds in physical units
lb = [gammaLower;
      shiftLower;
      globalBroadLower;
      phiLower];

ub = [gammaUpper;
      shiftUpper;
      globalBroadUpper;
      phiUpper];

% Parameter scaling for nonlinear optimization
thetaScale = [gammaScale;
              shiftScale;
              globalBroadScale;
              phiScale];

%% 5. Check our previous model is reproduced

r0 = profiledResidual(theta, basisFids, t, mask, sFit);

fprintf("\nInitial relative residual: %.2f%%\n", ...
    100*norm(r0));

% This should be close to the earlier 58.81% result.
if abs(norm(r0) - 0.5881) > 0.05
    warning("Initial fit differs from previous result. Check data conventions.");
end

%% 6. Staged nonlinear optimization

stageIndices = { ...
    M+(1:M), ...             % Stage 1: individual shifts
    [1:M, 2*M+1], ...        % Stage 2: linewidths + Gaussian
    1:(2*M+2) ...            % Stage 3: all nonlinear parameters
};

stageNames = [ ...
    "Individual frequency shifts", ...
    "Individual linewidths and Gaussian", ...
    "Full joint refinement"
];

stageHistory = struct([]);

for stage = 1:3

    fprintf("\n========== STAGE %d: %s ==========\n", ...
        stage, stageNames(stage));

    % Start Gaussian slightly positive when releasing it
    if stage == 2
        theta(2*M+1) = 3;
    end

    active = stageIndices{stage};

    % Fixed reference for parameters not active in this stage
    thetaBase = theta;

    % Optimize scaled variables
    x0 = theta(active) ./ thetaScale(active);
    xLB = lb(active) ./ thetaScale(active);
    xUB = ub(active) ./ thetaScale(active);

    objective = @(x) profiledResidual( ...
        replaceParameters( ...
            thetaBase, active, x .* thetaScale(active)), ...
        basisFids, t, mask, sFit);

    options = optimoptions("lsqnonlin", ...
        "Algorithm", "trust-region-reflective", ...
        "Display", "iter", ...
        "FiniteDifferenceType", "forward", ...
        "FiniteDifferenceStepSize", 1e-4, ...
        "TypicalX", ones(numel(active),1), ...
        "MaxIterations", 150, ...
        "MaxFunctionEvaluations", 10000, ...
        "FunctionTolerance", 1e-8, ...
        "StepTolerance", 1e-8, ...
        "OptimalityTolerance", 1e-6);

    [xHat, ~, ~, exitflag, solverOutput] = ...
        lsqnonlin(objective, x0, xLB, xUB, options);

    % Restore physical parameters
    theta(active) = xHat .* thetaScale(active);

    rStage = profiledResidual( ...
        theta, basisFids, t, mask, sFit);

    stageHistory(stage).name = stageNames(stage);
    stageHistory(stage).theta = theta;
    stageHistory(stage).relativeResidual = norm(rStage);
    stageHistory(stage).exitflag = exitflag;
    stageHistory(stage).solverOutput = solverOutput;

    fprintf("Stage %d relative residual: %.2f%%\n", ...
        stage, 100*norm(rStage));
end

%% 7. Extract final fitted parameters

gammaHat = theta(1:M);
deltaHat = theta(M+1:2*M);
gHat = theta(2*M+1);
phiHat = theta(2*M+2);

[rFinal, cHat, sHat, betaHat, baselineHat] = ...
    profiledResidual( ...
        theta, basisFids, t, mask, sFit);

relativeError = norm(rFinal);

fprintf("\n========== FINAL RESULT ==========\n");
fprintf("Relative residual: %.2f%%\n", 100*relativeError);
fprintf("Gaussian FWHM: %.3f Hz\n", gHat);
fprintf("Global phase: %.3f degrees\n", rad2deg(phiHat));

parameterTable = table( ...
    metaboliteNames(:), ...
    cHat(:), ...
    gammaHat(:), ...
    deltaHat(:), ...
    'VariableNames', { ...
        'Metabolite', 'Amplitude', ...
        'LorentzianHz', 'FrequencyShiftHz'});

disp(parameterTable);

%% 8. Save results in workspace

fitResult.theta = theta;
fitResult.concentrations = cHat;
fitResult.parameterTable = parameterTable;
fitResult.gaussianFWHM = gHat;
fitResult.phaseDegrees = rad2deg(phiHat);
fitResult.relativeResidual = relativeError;
fitResult.measured = sFit;
fitResult.fitted = sHat;
fitResult.ppm = ppmFit;
fitResult.residual = sFit - sHat;
fitResult.stageHistory = stageHistory;
fitResult.baselineCoefficients = betaHat;
fitResult.baseline = baselineHat;
fitResult.metaboliteFit = sHat - baselineHat;

%% 9. Plot fitted spectrum and residuals

figure("Name", "Full Nonlinear Spectral Fit");
tiledlayout(3,1);

nexttile;
plot(ppmFit, real(sFit), 'b'); hold on;
plot(ppmFit, real(sHat), 'r');
set(gca, 'XDir', 'reverse');
legend("Measured", "Nonlinear fit");
title("Real Spectrum");
grid on;

nexttile;
plot(ppmFit, imag(sFit), 'b'); hold on;
plot(ppmFit, imag(sHat), 'r');
set(gca, 'XDir', 'reverse');
legend("Measured", "Nonlinear fit");
title("Imaginary Spectrum");
grid on;

nexttile;
plot(ppmFit, real(sFit-sHat), 'b'); hold on;
plot(ppmFit, imag(sFit-sHat), 'r');
set(gca, 'XDir', 'reverse');
legend("Real residual", "Imaginary residual");
xlabel("Chemical shift (ppm)");
title(sprintf("Relative residual: %.2f%%", ...
    100*relativeError));
grid on;

sgtitle("P01 | Division 36 | Full Nonlinear Model");

figure;
plot(ppmFit, real(baselineHat)); hold on;
plot(ppmFit, imag(baselineHat));
set(gca, 'XDir', 'reverse');
legend('Real baseline', 'Imaginary baseline');
xlabel('Chemical shift (ppm)');
grid on;


%% LOCAL FUNCTIONS

function patients = getPatients(outputs, patientIDs)
    M = numel(patientIDs);
    patients = cell(1, M);
    allPatientIDs = string({outputs.patientResults.patientID});
    
    for id = 1:M
        p = find(allPatientIDs == patientIDs(id), 1);
        assert(~isempty(p), "Patient %s not found.", patientIDs(id));
        patients{id} = outputs.patientResults(p);
    end
end

function [fidSignals, specSignals, lengths, basisMatrices] = ...
    extractSignals(patientClasses)
    N = 10000;  % Pre-allocation maximum value, should later be trimmed using lengths
    M = numel(patientClasses);
    lengths = zeros(1, M);
    fidSignals = complex(zeros(N, M));
    specSignals = complex(zeros(N, M));
    basisMatrices = cell(1, M);
    
    for p = 1:M
        patient = patientClasses{p};
    
        T = patient.originalData.fileTable;
        idx = find(T.division == 36 & T.part == 1, 1);
    
        assert(~isempty(idx), "Requested measurement not found.");
    
        fid = patient.originalData.fids(:, idx);
        fidLen = numel(fid);

        lengths(p) = fidLen;
        fidSignals(1: fidLen, p) = fid;

        s = fftshift(fft(fid));
        specSignals(1:fidLen, p) = s; 

        basisMatrices{p} = getBasis(patient, fidLen);


    end
end

function B = getBasis(patient, n)
    basisStored = patient.preFitBasis.spectra;
    metaboliteNames = string(patient.preFitBasis.metaboliteNames);
    M = numel(metaboliteNames);
    
    B = complex(zeros(n, M));
    
    for j = 1:M
        if iscell(basisStored)
            bj = basisStored{j}(:);
        else
            bj = basisStored(:, j);
        end
    
        % Undo the original twofold zero filling
        bPadded = ifft(ifftshift(bj));
        assert(numel(bPadded) == 2*n);
    
        bFid = bPadded(1:n);
    
        % Transform to the measured frequency grid
        B(:, j) = fftshift(fft(bFid));
    end
end


function plotFIDSpec(outputs, patientID, metabName)
    patientToPlot = getPatients(outputs, patientID);
    [fid, S, length, basisMat] = extractSignals(patientToPlot);
    
    n = length(1);
    fid = fid(1:n);
    S = S(1:n);
    
    patientToPlot = patientToPlot{1};
    
    allNames = string(patientToPlot.preFitBasis.metaboliteNames);
    nameInd = find(allNames == metabName, 1);
    assert(~isempty(nameInd), metabName + " not found in basis.");
    bChosenMetab = basisMat{1}(:, nameInd);
    
    % Construct chemical shift axis
    hdr = patientToPlot.preFitBasis.globalHeader;
    
    dt = hdr.badelt;
    SW = 1 / dt;
    freqHz = -SW/2 + (0:n-1).' * (SW/n);
    
    B0 = hdr.hzpppm / ...
        SpinsJ.GetGyromagneticRatio("1H");
    
    ppm = SpinsJ.HzToPPM( ...
        -freqHz, B0, "1H", SpinsJ.waterPPM);
    
    mask = ppm >= 0.2 & ppm <= 4.0;
    
    % Plot everything
    
    figure("Name", "MRS Signal Inspection");
    tiledlayout(2,2);
    
    % Measured FID
    nexttile;
    t = (0:n-1).' * dt * 1000;
    plot(t, real(fid), 'b'); hold on;
    plot(t, imag(fid), 'r');
    xlabel("Time (ms)");
    ylabel("Signal");
    title("Measured Complex FID");
    legend("Real", "Imaginary");
    xlim([0 200]); grid on;
    
    % Measured complex spectrum
    nexttile;
    plot(ppm(mask), real(S(mask)), 'b'); hold on;
    plot(ppm(mask), imag(S(mask)), 'r');
    set(gca, 'XDir', 'reverse');
    xlabel("ppm");
    title("Measured Complex Spectrum");
    legend("Real", "Imaginary");
    grid on;
    
    % Measured magnitude vs pre-fit NAA basis
    nexttile;
    sNorm = abs(S) / max(abs(S(mask)));
    bNorm = abs(bChosenMetab) / max(abs(bChosenMetab(mask)));
    
    plot(ppm(mask), sNorm(mask), 'b'); hold on;
    plot(ppm(mask), bNorm(mask), 'r');
    set(gca, 'XDir', 'reverse');
    xlabel("ppm");
    title("Measured vs Pre-fit NAA (Magnitude)");
    legend("Measured", "NAA basis");
    grid on;
    
    % LCModel reference
    nexttile;
    T = patientToPlot.originalData.fileTable;
    idx = find(T.division == 36 & T.part == 1, 1);
    F = patientToPlot.fitData(idx);
    
    plot(F.axis, F.phasedData, 'b'); hold on;
    plot(F.axis, F.fitData, 'r');
    set(gca, 'XDir', 'reverse');
    xlabel("ppm");
    title("LCModel Result");
    legend("Measured", "LCModel fit");
    grid on;
    
    sgtitle("Patient P01 | Division 36 | Part 1");
end

function [sFit, BFit, ppmFit] = getAxis(patientClasses, spectrumSigs, ...
    sigLengths, basisMat)
    numPatients = numel(sigLengths);
    sFit = cell(1, numPatients);
    BFit = cell(1, numPatients);
    ppmFit = cell(1, numPatients);
        for patientID = 1: numPatients
            patient = patientClasses{patientID};
            n = sigLengths(patientID);
            s = spectrumSigs(1:n, patientID);
            B = basisMat{patientID};
            hdr = patient.preFitBasis.globalHeader;
        
            SW = 1 / hdr.badelt;
            freqHz = -SW/2 + (0:n-1).' * SW/n;
        
            B0 = hdr.hzpppm / SpinsJ.GetGyromagneticRatio("1H");
            ppm = SpinsJ.HzToPPM(-freqHz, B0, "1H", SpinsJ.waterPPM);
        
            % Select fitting range
            mask = ppm >= 0.2 & ppm <= 4.0;
        
            sFit{patientID} = s(mask);
            BFit{patientID} = B(mask, :);
            ppmFit{patientID} = ppm(mask);
        end
end

function [cHat, sHat, residual] = nnlsFit(BFit, sFit, ppmFit, opts)
    arguments
        BFit
        sFit
        ppmFit
        opts.Plot (1,1) logical = true
        opts.PatientIDs = []
    end
    
    numPatients = numel(BFit);
    sHat = cell(numPatients, 1);
    cHat = cell(numPatients, 1);
    residual = cell(numPatients, 1);
    
    if isempty(opts.PatientIDs)
        selectedPatients = 1:numel(BFit);
    else
        selectedPatients = opts.PatientIDs;
    end
    
    for p = selectedPatients
        % Extract patient data
        B = BFit{p};
        s = sFit{p};
        ppm = ppmFit{p};
    
        s = s(:);
        ppm = ppm(:);
    
        % Normalize basis columns for numerical stability    
        scale = sqrt(sum(abs(B).^2, 1));
        assert(all(isfinite(scale) & scale > 0), ...
            "Invalid or zero basis column.");
    
        Bnormalized = B ./ scale;
    
        % Stack real and imaginary components
        A = [real(Bnormalized); imag(Bnormalized)];
        y = [real(s); imag(s)];
    
        % Solve nonnegative least squares
        aHat = lsqnonneg(A, y);
    
        % Restore original basis scaling
        cHat{p} = aHat ./ scale(:);
    
        % Reconstruct fitted spectrum
        sHat{p} = B * cHat{p};
    
        % Calculate residual
        residual{p} = s - sHat{p};
    
        relativeError = norm(residual{p}) / norm(s);
    
        fprintf("Patient %d | Relative residual: %.2f%%\n", ...
            p, 100 * relativeError);
    
    
        % Optional plot measured, fitted, and residual spectra
        if opts.Plot
            figure("Name", "Linear Spectral Fit");
            tiledlayout(3,1);
    
            nexttile;
            plot(ppm, real(s), 'b'); hold on;
            plot(ppm, real(sHat{p}), 'r');
            set(gca, 'XDir', 'reverse');
            ylabel("Real");
            legend("Measured", "NNLS fit");
            title("Real Spectrum");
            grid on;
    
            nexttile;
            plot(ppm, imag(s), 'b'); hold on;
            plot(ppm, imag(sHat{p}), 'r');
            set(gca, 'XDir', 'reverse');
            ylabel("Imaginary");
            title("Imaginary Spectrum");
            grid on;
    
            nexttile;
            plot(ppm, real(residual{p}), 'b'); hold on;
            plot(ppm, imag(residual{p}), 'r');
            set(gca, 'XDir', 'reverse');
            xlabel("Chemical shift (ppm)");
            ylabel("Residual");
            legend("Real", "Imaginary");
            title(sprintf("Residual | Relative error: %.1f%%", ...
                100*relativeError));
            grid on;
    
            sgtitle(sprintf("Patient %d | Initial Linear Fit", p));
        end
    end
end

function theta = replaceParameters(theta, active, values)
    theta(active) = values(:);
end

function [r, cHat, sHat, betaHat, baselineHat] = profiledResidual( ...
    theta, basisFids, t, mask, sFit)

    M = size(basisFids, 2);

    % Extract nonlinear parameters
    gamma = theta(1:M);
    delta = theta(M+1:2*M);
    g = theta(2*M+1);
    phi = theta(2*M+2);

    % Individual Lorentzian broadening and frequency shifts
    lorentzShift = exp( ...
        (-pi*t)*gamma.' + ...
        (1i*2*pi*t)*delta.');

    % Shared Gaussian broadening
    gaussian = exp( ...
        -((pi*g*t).^2)/(4*log(2)));

    % Modify the original time-domain basis
    modifiedFids = basisFids .* lorentzShift .* gaussian;

    % FFT and global phase
    Bfull = exp(1i*phi) * ...
        fftshift(fft(modifiedFids, [], 1), 1);

    Bfit = Bfull(mask, :);

    % Two fixed-shape baseline functions
    K = numel(sFit);
    p = linspace(-1, 1, K).';

    baselineBasis = [ ...
        ones(K,1), ...
        p, ...
        1i*ones(K,1), ...
        1i*p];

    % Normalize basis columns
    scales = sqrt(sum(abs(Bfit).^2, 1));

    assert(all(isfinite(scales) & scales > 0), ...
        "Invalid basis columns.");

    Bnorm = Bfit ./ scales;

    % Complex NNLS converted to a real problem
    
    % Stack real and imaginary components
    A = [real(Bnorm); imag(Bnorm)];
    y = [real(sFit); imag(sFit)];
    
    % Stack the two baseline columns
    D = [real(baselineBasis); imag(baselineBasis)];
    
    % Remove components that the baseline can explain
    Aprojected = A - D * (D \ A);
    yprojected = y - D * (D \ y);
    
    % Estimate nonnegative metabolite amplitudes
    aHat = lsqnonneg(Aprojected, yprojected);
    
    % Estimate signed baseline coefficients
    betaHat = D \ (y - A*aHat);
    
    % Restore metabolite amplitudes
    cHat = aHat ./ scales(:);
    
    % Reconstruct both contributions
    metaboliteFit = Bfit * cHat;
    baselineHat = baselineBasis * betaHat;
    
    sHat = metaboliteFit + baselineHat;


    % Normalized real residual for lsqnonlin
    difference = sFit - sHat;

    r = [real(difference); imag(difference)] / norm(sFit);
end
