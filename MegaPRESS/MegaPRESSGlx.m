clear all;
close all;
clc;

%% Configuration
basisDir = "C:\Users\User\Documents\Thesis_Lab\MRS_Project\MegaPRESS";

% Each row defines [lowerPPM, upperPPM]
ppmRanges = [1.5, 4.2];

% Assumed complex noise variance per selected frequency bin
sigma2 = 1;

%% Load MEGA-PRESS basis functions
[BONAll, headersON, ppmON, ~] = ...
    VDILCM.ReadBasis(fullfile(basisDir, "editON.basis"));

[BOFFAll, headersOFF, ppmOFF, ~] = ...
    VDILCM.ReadBasis(fullfile(basisDir, "editOFF.basis"));

whos BON BOFF ppmON ppmOFF

disp(class(BONAll))
disp(headersON)

%% Construct Glu/Gln basis matrices

% Columns: [Glu, Gln]
B_ON  = [BONAll{1}(:),  BONAll{2}(:)];
B_OFF = [BOFFAll{1}(:), BOFFAll{2}(:)];

% Difference basis
B_DIFF = B_ON - B_OFF;

% Shared ppm axis
ppm = ppmON(:);

% Verify dimensions
assert(isequal(size(B_ON), size(B_OFF)));
assert(size(B_ON, 1) == numel(ppm));

disp(size(B_ON))

%% Select spectral regions of interest

% Construct a logical mask
roiMask = false(size(ppm));

for k = 1:size(ppmRanges, 1)
    roiMask = roiMask | ...
        (ppm >= ppmRanges(k,1) & ...
        ppm <= ppmRanges(k,2));
end

% Restrict the basis matrices
B_ON_ROI   = B_ON(roiMask, :);
B_OFF_ROI  = B_OFF(roiMask, :);
B_DIFF_ROI = B_DIFF(roiMask, :);

ppmROI = ppm(roiMask);

fprintf('Selected %d out of %d spectral points.\n', ...
    nnz(roiMask), numel(ppm));

%% CRLB - OFF acquisition

% Select native frequency bins within the ROI
nativeRows = (1:2:numel(ppm))';
selectedRows = nativeRows(roiMask(nativeRows));

J_OFF = B_OFF(selectedRows, :);

% Fisher information
F_OFF = statisticalTools.fisherInformation(J_OFF, sigma2);

% CRLB covariance matrix
[C_OFF, infoOFF] = statisticalTools.crlb(F_OFF);

% Display results
disp('Fisher Information:')
disp(F_OFF)

disp('CRLB covariance:')
disp(C_OFF)

disp('CRLB standard deviations [Glu, Gln]:')
disp(sqrt(diag(C_OFF)))

%% CRLB - ON acquisition

% Jacobian
J_ON = B_ON(selectedRows, :);

% Fisher information
F_ON = statisticalTools.fisherInformation(J_ON, sigma2);

% CRLB covariance
[C_ON, infoON] = statisticalTools.crlb(F_ON);

% Compare standard deviation bounds
sdOFF = sqrt(diag(C_OFF));
sdON  = sqrt(diag(C_ON));

results = table( ...
    ["Glu"; "Gln"], sdOFF, sdON, sdON ./ sdOFF, ...
    'VariableNames', {'Metabolite', 'OFF_SD', 'ON_SD', 'ON_to_OFF'});

disp(results)

%% CRLB - Joint ON + OFF

J_JOINT = [J_ON; J_OFF];

F_JOINT = statisticalTools.fisherInformation(J_JOINT, sigma2);

[C_JOINT, infoJOINT] = statisticalTools.crlb(F_JOINT);

sdJOINT = sqrt(diag(C_JOINT));

results = table( ...
    ["Glu"; "Gln"], sdOFF, sdON, sdJOINT, ...
    'VariableNames', {'Metabolite', 'OFF_SD', 'ON_SD', 'JOINT_SD'});

disp(results)

%% CRLB - Joint with equal acquisition time

F_JOINT_equal = F_JOINT / 2;

[C_JOINT_equal, ~] = statisticalTools.crlb(F_JOINT_equal);

sdJOINT_equal = sqrt(diag(C_JOINT_equal));

results = table( ...
    ["Glu"; "Gln"], sdOFF, sdJOINT_equal, ...
    sdJOINT_equal ./ sdOFF, ...
    'VariableNames', ...
    {'Metabolite', 'OFF_SD', 'JOINT_Equal_SD', 'Ratio'});

disp(results)

%% CRLB - Glx concentration sum

w = [1; 1];

% CRLB variance for Glx
varGlxOFF = w' * C_OFF * w;
varGlxJOINT = w' * C_JOINT_equal * w;

% Standard deviation bounds
sdGlxOFF = sqrt(varGlxOFF);
sdGlxJOINT = sqrt(varGlxJOINT);

fprintf('Glx OFF SD: %.6f\n', sdGlxOFF);
fprintf('Glx JOINT SD: %.6f\n', sdGlxJOINT);
fprintf('JOINT / OFF: %.4f\n', sdGlxJOINT / sdGlxOFF);
