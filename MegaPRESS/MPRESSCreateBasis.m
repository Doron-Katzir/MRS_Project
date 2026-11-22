clear all
close all
clc

%% Define spin, sequence and post-processing variables
% Spin system parameters
B0 = 3; % Tesla
concGlu = 1; % mM.
concGln = 1;
concCr = 1;

% Sequence parameters
% The typical MEGAPRESS echo time is 68 ms. This ensures that the GABA
% peak at 3 ppm yields a nice looking "pseudo-doublet". Try changing it 
% and see how it affects the results!
TE = 68; % ms. 

% Plotting parameters
FWHMHz = 5; % Hz

% Parameters for exporting the basis set
basisFWHM = 1; % Hz. Keep it narrow for fitting.
dwellTime = 1e-3; % Sec
numAcqPts = 512; 

% Define metabolites
metabolites = [
    "Glu", "Gln", "Pcr"
    ];

spins = SpinsJ('B0', B0);

for m = 1:numel(metabolites)
    spins.AddMetab(metabolites(m), ...
        "concentration", 1);
end

%% Create the MEGAPRESS sequence
seqON  = Sequence.GetIdealSequence("MEGAPRESSON", "TE", TE);
seqOFF = Sequence.GetIdealSequence("MEGAPRESSOFF", "TE", TE);

%% Acquire data
[spinsOutON,  TTON]  = seqON.Apply (spins, "isVerbose", true);
[spinsOutOFF, TTOFF] = seqOFF.Apply(spins, "isVerbose", true);

TTDIFF = TTON - TTOFF;

%% Plotting the resulting spectra.
TTON.Plot  ("FWHM", FWHMHz, "B0", B0);
TTOFF.Plot ("FWHM", FWHMHz, "B0", B0);
TTDIFF.Plot("FWHM", FWHMHz, "B0", B0);

%% Create and export basis functions

outputDir = "C:\Users\User\Documents\Thesis_Lab\MRS_Project\MegaPRESS";

if ~isfolder(outputDir)
    mkdir(outputDir);
end

basisON = VDIBasis(TTON, ...
    "dwellTime", dwellTime, "B0", B0, "numAcqPts", numAcqPts);

basisOFF = VDIBasis(TTOFF, ...
    "dwellTime", dwellTime, "B0", B0, "numAcqPts", numAcqPts);

basisDIFF = VDIBasis(TTDIFF, ...
    "dwellTime", dwellTime, "B0", B0, "numAcqPts", numAcqPts);

basisON.ExportBasisToLCModel( ...
    fullfile(outputDir, "editON.basis"), ...
    "TE", TE, "FWHM", basisFWHM, "ampDSS", 0);

basisOFF.ExportBasisToLCModel( ...
    fullfile(outputDir, "editOFF.basis"), ...
    "TE", TE, "FWHM", basisFWHM, "ampDSS", 0);

basisDIFF.ExportBasisToLCModel( ...
    fullfile(outputDir, "editDIFF.basis"), ...
    "TE", TE, "FWHM", basisFWHM, "ampDSS", 0);

