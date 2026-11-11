cfg = ProjectConfig();

cfg.paths.rootDir = "C:\Users\User\Documents\Thesis_Lab";
cfg.paths.dataDir = fullfile(cfg.paths.rootDir, "Data");
cfg.paths.coordDir = fullfile(cfg.paths.rootDir, "LCMFit");

cfg.input.mode = "singleFile";
cfg.input.singleFile = fullfile( ...
    cfg.paths.dataDir, ...
    "meas_MID00020_FID54986_eja_svs_slaser_TE_80_1_PRE.dat");

cfg.input.singlePatientID = "P01";
cfg.output.deleteOldLCModelFiles = false;

processedPatients = Splice_data_multi_patient_input(cfg);