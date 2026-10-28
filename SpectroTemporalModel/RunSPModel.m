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
