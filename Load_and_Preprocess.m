data_file = "meas_MID00090_FID32072_eja_svs_slaser_TE_80_r0.dat";
data_path = "C:\Users\User\Documents\Thesis_Lab\Data\" + data_file;
basis_file = "sLaser_TE70_20x20x20.basis";
basis_path = "C:\Users\User\Documents\Thesis_Lab\Data\" + basis_file;

%%
data = VDIIO.LoadTwix(data_path, "isICEChop", true);

disp(data.table);
data.AddCoils;

num_points = 2 ^ floor(log2(data.numSpecPts));
data.ChopPts("numPts", num_points);
disp(data.table);

data.FT;
disp(data.table);

data.Average("dim", "set");
disp(data.table);

data.Plot;

VDIPlot.PlotLCMBasis(basis_path);

%%
data.FitLCModel(basis_path, 'isVerbose', true);
VDIPlot.Plot