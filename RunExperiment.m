clear;
clc;

cfg = ProjectConfig();

%% User input for this run

cfg.paths.rootDir = "C:\Users\doronkatzir1\Desktop\Thesis_Lab";
cfg.paths.dataDir = fullfile(cfg.paths.rootDir, "Data");
cfg.paths.coordDir = fullfile(cfg.paths.rootDir, "LCMFit");

%% Splice data and save coord files
% Splice only one file

% cfg.input.mode = "singleFile";
% cfg.input.singlePatientID = "P11";
% cfg.input.singleFileName = "meas_MID00090_FID32072_eja_svs_slaser_TE_80_r0.dat";
% cfg.input.singleFile = fullfile(cfg.paths.dataDir, cfg.input.singleFileName);
% processedPatients = Splice_data_multi_patient_input(cfg);

% Splice all files in the working directory

% cfg.input.mode = "directory";
% 
% cfg.input.directory = cfg.paths.dataDir;
% cfg.input.filePattern = "*.dat";
% cfg.input.recursive = false;
% 
% processedPatients = Splice_data_multi_patient_input(cfg);

%% Load coord files and analyze
% Load only one fitted patient folder

% cfg.load.mode = "singleSubfolder";
% cfg.load.selectedPatientID = "P11";
% cfg.load.coordDir = cfg.paths.coordDir;
% 
% outputs = Load_subsets_multi_patient_input(cfg);

% Load all subfolders

cfg.load.mode = "allSubfolders";
cfg.load.coordDir = cfg.paths.coordDir;

outputs = Load_subsets_multi_patient_input(cfg);