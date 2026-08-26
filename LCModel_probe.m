% Load basis file
basisFile = "C:\Users\User\Documents\Thesis_Lab\LCMFit\P01\For_division_P01_meas_MID00020_FID54986_eja_svs_slaser_TE_80_1_PRE_1_.basis";
[basisCell, basisHdr, ppmBasis, globalHdr] = VDILCM.ReadBasis(basisFile);
% basisCell contains the actual pre-fit basis functions
% basisHdr is metadata on the basis
% ppmBasis is the x axis

% Plot one metab pre-fit basis function
metabInd = 1;
plot(ppmBasis, real(basisCell{metabInd}));
set(gca,'XDir','reverse');
title("Basis function for " + basisHdr(1).metabo);
xlabel("ppm");
ylabel("Real Amplitude");

% Load coord file
coordFile = "C:\Users\User\Documents\Thesis_Lab\LCMFit\P01\P01_Division_1_part_1.basis.coord";
[quant, fitData] = VDIIO.ReadLCMCoord(coordFile);

% Plot one metab contribution (including scaling)
j = find(fitData.basisMetName == "NAA");

figure;
naaPostFit = fitData.basisData(:,j);
plot(fitData.axis, naaPostFit)
set(gca,'XDir','reverse')
xlabel('ppm')

% Plot one metab post-fit basis
naaContribution = fitData.basisData(:,j);
cNAA = 1.06e-7;

naaPostFitCandidate = naaPostFit / cNAA;
plot(fitData.axis, naaPostFitCandidate)
set(gca,'XDir','reverse')
xlabel('ppm')

% Load print file
printFile = "C:\Users\User\Documents\Thesis_Lab\LCMFit\P01\P01_Division_1_part_1.basis.print";