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

% 1. Select patient and averaged measurement

patientIDs = string({outputs.patientResults.patientID});
p = find(patientIDs == "P01", 1);

assert(~isempty(p), "Patient P01 not found.");
patient = outputs.patientResults(p);

T = patient.originalData.fileTable;
idx = find(T.division == 36 & T.part == 1, 1);

assert(~isempty(idx), "Requested measurement not found.");

fid = patient.originalData.fids(:, idx);
n = numel(fid);

% 2. Convert measured FID to spectrum

S = fftshift(fft(fid));

% 3. Load original NAA basis and undo zero filling

names = string(patient.preFitBasis.metaboliteNames);
j = find(names == "NAA", 1);

assert(~isempty(j), "NAA not found in basis.");

basisStored = patient.preFitBasis.spectra;

if iscell(basisStored)
    b8192 = basisStored{j}(:);
else
    b8192 = basisStored(:, j);
end

% Reverse the export-time zero filling
bFidPadded = ifft(ifftshift(b8192));
assert(numel(bFidPadded) == 2*n);

bFid = bFidPadded(1:n);
Bnaa = fftshift(fft(bFid));

% Check removed tail is approximately zero
tailRatio = norm(bFidPadded(n+1:end)) / ...
            max(norm(bFid), eps);

fprintf("Zero-filled tail ratio: %.3e\n", tailRatio);

% 4. Construct chemical shift axis

hdr = patient.preFitBasis.globalHeader;

dt = hdr.badelt;
SW = 1 / dt;

freqHz = -SW/2 + (0:n-1).' * (SW/n);

B0 = hdr.hzpppm / ...
    SpinsJ.GetGyromagneticRatio("1H");

ppm = SpinsJ.HzToPPM( ...
    -freqHz, B0, "1H", SpinsJ.waterPPM);

mask = ppm >= 0.2 & ppm <= 4.0;

% 5. Plot everything

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
bNorm = abs(Bnaa) / max(abs(Bnaa(mask)));

plot(ppm(mask), sNorm(mask), 'b'); hold on;
plot(ppm(mask), bNorm(mask), 'r');
set(gca, 'XDir', 'reverse');
xlabel("ppm");
title("Measured vs Pre-fit NAA (Magnitude)");
legend("Measured", "NAA basis");
grid on;

% LCModel reference
nexttile;
F = patient.fitData(idx);

plot(F.axis, F.phasedData, 'b'); hold on;
plot(F.axis, F.fitData, 'r');
set(gca, 'XDir', 'reverse');
xlabel("ppm");
title("LCModel Result");
legend("Measured", "LCModel fit");
grid on;

sgtitle("Patient P01 | Division 36 | Part 1");


%% Build the complete pre-fit basis matrix
s = fftshift(fft(fid));
basisStored = patient.preFitBasis.spectra;
metaboliteNames = string(patient.preFitBasis.metaboliteNames);
M = numel(metaboliteNames);

B = zeros(n, M);

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

% Construct ppm axis from acquisition metadata
hdr = patient.preFitBasis.globalHeader;

SW = 1 / hdr.badelt;
freqHz = -SW/2 + (0:n-1).' * SW/n;

B0 = hdr.hzpppm / SpinsJ.GetGyromagneticRatio("1H");
ppm = SpinsJ.HzToPPM(-freqHz, B0, "1H", SpinsJ.waterPPM);

% Select fitting range
mask = ppm >= 0.2 & ppm <= 4.0;

sFit = s(mask);
BFit = B(mask, :);
ppmFit = ppm(mask);

fprintf("Measured spectrum: %d points\n", numel(sFit));
fprintf("Basis matrix: %d x %d\n", size(BFit));


%% Complex nonnegative least squares

% Normalize basis columns for numerical stability
scale = sqrt(sum(abs(BFit).^2, 1));
assert(all(isfinite(scale) & scale > 0), ...
    "Invalid or zero basis column.");

Bnormalized = BFit ./ scale;

% Stack real and imaginary components
A = [real(Bnormalized); imag(Bnormalized)];
y = [real(sFit); imag(sFit)];

% Solve nonnegative least squares
aHat = lsqnonneg(A, y);

% Restore original basis scaling
cHat = aHat ./ scale(:);

% Reconstruct fitted spectrum
sHat = BFit * cHat;

% Calculate residual
residual = sFit - sHat;

relativeError = norm(residual) / norm(sFit);

fprintf("Relative residual: %.2f%%\n", ...
    100 * relativeError);



%% Plot measured, fitted, and residual spectra

figure("Name", "First Linear Spectral Fit");
tiledlayout(3,1);

nexttile;
plot(ppmFit, real(sFit), 'b'); hold on;
plot(ppmFit, real(sHat), 'r');
set(gca, 'XDir', 'reverse');
ylabel("Real");
legend("Measured", "NNLS fit");
title("Real Spectrum");
grid on;

nexttile;
plot(ppmFit, imag(sFit), 'b'); hold on;
plot(ppmFit, imag(sHat), 'r');
set(gca, 'XDir', 'reverse');
ylabel("Imaginary");
title("Imaginary Spectrum");
grid on;

nexttile;
plot(ppmFit, real(residual), 'b'); hold on;
plot(ppmFit, imag(residual), 'r');
set(gca, 'XDir', 'reverse');
xlabel("Chemical shift (ppm)");
ylabel("Residual");
legend("Real", "Imaginary");
title(sprintf("Residual | Relative error: %.1f%%", ...
    100*relativeError));
grid on;

sgtitle("P01 | Division 36 | Initial Linear Fit");



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
