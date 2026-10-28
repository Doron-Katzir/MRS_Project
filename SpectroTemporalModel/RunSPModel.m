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

%% 1. Select patient and averaged measurement

patientIDs = string({outputs.patientResults.patientID});
p = find(patientIDs == "P01", 1);

assert(~isempty(p), "Patient P01 not found.");
patient = outputs.patientResults(p);

T = patient.originalData.fileTable;
idx = find(T.division == 36 & T.part == 1, 1);

assert(~isempty(idx), "Requested measurement not found.");

fid = patient.originalData.fids(:, idx);
n = numel(fid);

%% 2. Convert measured FID to spectrum

S = fftshift(fft(fid));

%% 3. Load original NAA basis and undo zero filling

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

%% 4. Construct chemical shift axis

hdr = patient.preFitBasis.globalHeader;

dt = hdr.badelt;
SW = 1 / dt;

freqHz = -SW/2 + (0:n-1).' * (SW/n);

B0 = hdr.hzpppm / ...
    SpinsJ.GetGyromagneticRatio("1H");

ppm = SpinsJ.HzToPPM( ...
    -freqHz, B0, "1H", SpinsJ.waterPPM);

mask = ppm >= 0.2 & ppm <= 4.0;

%% 5. Plot everything

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
