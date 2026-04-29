img = VDIIO.LoadTwix('meas_MID00090_FID32072_eja_svs_slaser_TE_80_r0.dat', 'isICEChop', true);
img.AddCoils;
img.ChopPts("numPts", 4096);
img.FT;
img.Phase; % This will automatically (hopefully!) phase your data
% img.Average("dim", "set");

%% Creating a simple basis set
% We will simulate a .basis file with a set of basis functions (metabolites)
% using the VDI routines with some explanations along the way.
% Simplest case: "single molecule" at x=y=z=0, with an idealized simple
% sequence (usually PRESS)
% All of the simulations are carried out using the SpinsJ class.
spins = SpinsJ('B0', img.B0);
% Let's populate 18 different metabolite types, which will be stored in
% the spins.metab property
metabList = ["NAA", "NAAG", "Cr", "PCr", "GPC", "PCh", "Glu", "Gln", "GABA", "GSH", "Tau", "Asc", "Glc", "Ace", "mI", "sI", "Asp", "Lac"];
spins.AddMetab(metabList);
% Now will create an idealized PRESS sequence.
% Use the actual TE from the header data in your VDIImageND object,
% which can be read using the GetTwixHeaderReport method.
TE = VDIIO.GetTwixHeaderReport(img.metadata.hdr).TE;
seq = Sequence.GetIdealSequence("PRESS", "TE", TE);
[spinsOut, TT] = seq.Apply(spins, "isVerbose", true);
% The output will be a 1x18 array of TransitionTable objects, each of
% which will contain the "results" of the simulation for each of the 18
% metabolites. For example, TT(4) is for PCr:
% freqHz amp phaseDeg PPM3T PPM7T
% _______ ____ __________ _______ ______
%
% -490.67 0.75 -0.011245 0.83862 3.0337
% -222.9 0.5 -0.0051084 2.935 3.9321
% Create a VDI basis-set object
myBasis = VDIBasis(TT, 'B0', img.B0, ...
    'numAcqPts', img.numSpecPts, ...
    'dwellTime', img.dwellTime);
myBasis.ExportBasisToLCModel("Doron.basis", 'TE', TE);

%% Subdividing a dataset
% You will need to manually subdivide img.data and create individual
% VDIImageND objects, which will then be fit (all using the same basis set).
img1 = img.Copy;
img1.data = img1.data(:,:,:,:,1:18);
img1.FitLCModel("Doron.basis", 'isVerbose', true);
% TODO: Look at time courses of metabolites at different "temporal resolutions"
% 1x36
% 2x18
% 4x9
% 6x6
% 9x4
% 18x2
% 36x1
% For each you can plot the same time course (over the same x-axis time range!)