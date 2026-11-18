%% ========================================================================
%             MEGAPRESS Simulations: From sequence to spectrum
% =========================================================================
%
% In this example we will simulate a MEGAPRESS experiment, with ON and OFF
% acquisitions, acquire spectra, and demonstrate how a difference spectrum
% should look like.
%
% =========================================================================

clear all
close all
clc

%% ========================================================================
% Define spin, sequence and post-processing variables
% =========================================================================

% Spin system parameters
B0 = 3; % Tesla
concGABA = 1; % mM. Typical in-vivo concentration.
concCr = 5; % mM. Typical in-vivo concentration.

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


%% ========================================================================
% Create the spins. 
%
% For this simulation we will initialize a sample with one mM of GABA and
% 5 mM of creatine. This will enable us to observe how GABA and Cr fall
% on top of each other in a simple acquisition, and how editing eliminates
% the non-coupled creatine peaks.
% =========================================================================

spins = SpinsJ('B0', B0);
spins.AddMetab("GABA",  "concentration", concGABA);
spins.AddMetab("Cr", "concentration", concCr);

%% ========================================================================
% Create the MEGAPRESS sequence
%
% There are actually two sequences: "ON" and "OFF", which differ by the
% frequency at which editing pulses are applied. In this simplistic 
% implementation:
% 1. Both editing pulses will be modeled as hard pulses which simply
%    act on the [1.7, 2.1] ppm range, affected the GABA resonance at 1.9 ppm
%    and thus indirectly affecting the two other GABA resonances.
% 2. In the OFF acquisition the pulses will be completely removed (this is
%    unlike a "real" sequence where the frequency range is simply shifted
%    to some region of the spectrum where there is "nothing", e.g. 8-9 ppm).
%
%                          ON/OFF                      ON/OFF
%      90       180        180(S)         180          180(S)
%       _        _           _             _             _         
%      | |      | |         | |           | |           | |       
%      | |      | |         | |           | |           | |         
% RF __| |______| |_________| |___________| |___________| |________Acquire
%
% =========================================================================

seqON  = Sequence.GetIdealSequence("MEGAPRESSON", "TE", TE);
seqOFF = Sequence.GetIdealSequence("MEGAPRESSOFF", "TE", TE);

% Note: you can type 
% >> seqON.Plot;
% to see the sequence. You can also simply type "seqON" at the command
% line to see a summary table of all pulses.

%% ========================================================================
% Acquire data
% 
% The structures of interest here are TTON and TTOFF, which are the
% "transition tables" of the acquisition. Namely, they contain information
% about each observed line: its frequency, amplitude, and phase. If you
% type "TTON" at the command line you will get, e.g., a collection of all
% such lines (for each metabolite).
% =========================================================================

[spinsOutON,  TTON]  = seqON.Apply (spins, "isVerbose", true);
[spinsOutOFF, TTOFF] = seqOFF.Apply(spins, "isVerbose", true);

%% ========================================================================
% Post processing: create the difference spectrum (DIFF)
% =========================================================================

TTDIFF = TTON - TTOFF;

%% ========================================================================
% Plotting the resulting spectra.
% =========================================================================

TTON.Plot  ("FWHM", FWHMHz, "B0", B0);
TTOFF.Plot ("FWHM", FWHMHz, "B0", B0);
TTDIFF.Plot("FWHM", FWHMHz, "B0", B0);

%% ========================================================================
% Create and export .basis files
% =========================================================================

basisON   = VDIBasis(TTON,   "dwellTime", dwellTime, "B0", B0, "numAcqPts", numAcqPts);
basisOFF  = VDIBasis(TTOFF,  "dwellTime", dwellTime, "B0", B0, "numAcqPts", numAcqPts);
basisDIFF = VDIBasis(TTDIFF, "dwellTime", dwellTime, "B0", B0, "numAcqPts", numAcqPts);

basisON.ExportBasisToLCModel  ("editON.basis",   "TE", TE, "FWHM", basisFWHM, "ampDSS", 0);
basisOFF.ExportBasisToLCModel ("editOFF.basis",  "TE", TE, "FWHM", basisFWHM, "ampDSS", 0);
basisDIFF.ExportBasisToLCModel("editDIFF.basis", "TE", TE, "FWHM", basisFWHM, "ampDSS", 0);