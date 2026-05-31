%% Project Configuration Parameters

function cfg = ProjectConfig()
    rootDir = "C:\Users\User\Documents\Thesis_Lab";
    
    cfg.paths.rootDir = rootDir;
    cfg.paths.dataDir = fullfile(rootDir, "Data");
    cfg.paths.coordDir = fullfile(rootDir, "LCMFit");
    
    cfg.files.twixFile = fullfile(cfg.paths.dataDir, ...
        "meas_MID00090_FID32072_eja_svs_slaser_TE_80_r0.dat");
    
    cfg.preprocessing.isICEChop = true;
    cfg.preprocessing.numPts = 4096;
    
    cfg.subsets.setSizes = [1, 2, 4, 6, 9, 12, 18, 36];
    
    cfg.metabolites.basis = ["NAA", "NAAG", "Cr", "PCr", "GPC", "PCh", ...
        "Glu", "Gln", "GABA", "GSH", "Tau", "Asc", "Glc", "Ace", ...
        "mI", "sI", "Asp", "Lac"];
    
    cfg.metabolites.sum = ["GPC+PCh", "NAA+NAAG", "Cr+PCr", "Glu+Gln"];
    
    cfg.metabolites.analysis = [ ...
        cfg.metabolites.basis, ...
        cfg.metabolites.sum ...
        ];
    
    cfg.lcmodel.valueColumn = "sig";
    cfg.lcmodel.outputNamePattern = "Division_%d_part_%d.basis";
end