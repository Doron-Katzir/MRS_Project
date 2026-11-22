clear all;
close all;
clc;

%% Configuration
basisDir = "C:\Users\User\Documents\Thesis_Lab\MRS_Project\MegaPRESS";

% Each row defines [lowerPPM, upperPPM]
ppmRanges = [1.8, 2.6];

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
metNames = string({headersON.metabo});

assert(isequal(metNames, string({headersOFF.metabo})), ...
    'ON and OFF metabolite names do not match.');

% Load all basis functions, including PCr
B_ON = cell2mat(cellfun(@(x) x(:), BONAll, ...
    'UniformOutput', false));

B_OFF = cell2mat(cellfun(@(x) x(:), BOFFAll, ...
    'UniformOutput', false));

% Find target metabolite indices
idxGlu = find(metNames == "Glu");
idxGln = find(metNames == "Gln");

assert(isscalar(idxGlu) && isscalar(idxGln), ...
    'Glu or Gln missing from the basis set.');

% Shared ppm axis
ppm = ppmON(:);

% Verify dimensions before selecting metabolites
assert(isequal(size(B_ON), size(B_OFF)));
assert(size(B_ON, 1) == numel(ppm));

% Restrict CRLB analysis to Glu and Gln

idxKeep = [idxGlu, idxGln];

B_ON = B_ON(:, idxKeep);
B_OFF = B_OFF(:, idxKeep);

% Update names and indices
metNames = metNames(idxKeep);

idxGlu = 1;
idxGln = 2;

% Difference basis
B_DIFF = B_ON - B_OFF;

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

%% CRLB - ON acquisition

J_ON = B_ON(selectedRows, :);

F_ON = statisticalTools.fisherInformation(J_ON, sigma2);

[C_ON, infoON] = statisticalTools.crlb(F_ON);

%% CRLB - Joint ON + OFF

J_JOINT = [J_ON; J_OFF];

F_JOINT = statisticalTools.fisherInformation(J_JOINT, sigma2);

[C_JOINT, infoJOINT] = statisticalTools.crlb(F_JOINT);

%% CRLB - Joint with equal acquisition time

F_JOINT_equal = F_JOINT / 2;

[C_JOINT_equal, ~] = statisticalTools.crlb(F_JOINT_equal);

sdJOINT_equal = sqrt(diag(C_JOINT_equal));

%% CRLB - OFF + DIFF with equal acquisition time

% Difference spectrum: ON - OFF
J_DIFF = J_ON - J_OFF;

% Concatenated Jacobian
J_OFF_DIFF = [J_OFF; J_DIFF];

% Number of frequency samples
K = size(J_OFF, 1);
I = eye(K);

% Complex noise covariance matrix
Sigma_OFF_DIFF = 2 * sigma2 * [
    I, -I;
    -I, 2*I
    ];

% Fisher information
F_OFF_DIFF = statisticalTools.fisherInformation( ...
    J_OFF_DIFF, Sigma_OFF_DIFF);

% CRLB covariance matrix
[C_OFF_DIFF, infoOFF_DIFF] = ...
    statisticalTools.crlb(F_OFF_DIFF);

%% Display Fisher Information and CRLB matrices

% Metabolite labels
labels = cellstr(metNames);

% All acquisition strategies
caseNames = {
    'OFF'
    'ON'
    'JOINT'
    'JOINT - Equal Time'
    'OFF + DIFF - Equal Time'
};

F_matrices = {
    F_OFF
    F_ON
    F_JOINT
    F_JOINT_equal
    F_OFF_DIFF
};

C_matrices = {
    C_OFF
    C_ON
    C_JOINT
    C_JOINT_equal
    C_OFF_DIFF
};

% Display labeled matrices
for k = 1:numel(caseNames)

    fprintf('\n===== %s =====\n', caseNames{k});

    % Fisher Information Matrix
    F_table = array2table(F_matrices{k}, ...
        'VariableNames', labels, ...
        'RowNames', labels);

    disp('Fisher Information Matrix:')
    disp(F_table)

    % CRLB Covariance Matrix
    C_table = array2table(C_matrices{k}, ...
        'VariableNames', labels, ...
        'RowNames', labels);

    disp('CRLB Covariance Matrix:')
    disp(C_table)

end


%% Results - Glu, Gln and Glx
% Target metabolite indices
idxTarget = [idxGlu; idxGln];

% Standard deviation lower bounds for all metabolites
sdOFF = sqrt(diag(C_OFF));
sdON = sqrt(diag(C_ON));
sdJOINT = sqrt(diag(C_JOINT));
sdJOINT_equal = sqrt(diag(C_JOINT_equal));
sdOFF_DIFF = sqrt(diag(C_OFF_DIFF));

% Extract Glu and Gln results
results = table( ...
    metNames(idxTarget)', ...
    sdOFF(idxTarget), ...
    sdON(idxTarget), ...
    sdJOINT_equal(idxTarget), ...
    sdOFF_DIFF(idxTarget), ...
    sdJOINT_equal(idxTarget) ./ sdOFF(idxTarget), ...
    sdOFF_DIFF(idxTarget) ./ sdOFF(idxTarget), ...
    'VariableNames', ...
    {'Metabolite', 'OFF_SD', 'ON_SD', ...
     'JOINT_Equal_SD', 'OFF_DIFF_SD', ...
     'JOINT_Ratio', 'OFF_DIFF_Ratio'});

disp(results)

%% Glx CRLB

% Linear combination: Glx = Glu + Gln
w = zeros(numel(metNames), 1);

w(idxGlu) = 1;
w(idxGln) = 1;

% Variance lower bounds
varGlxOFF = w' * C_OFF * w;
varGlxJOINT = w' * C_JOINT_equal * w;
varGlxOFF_DIFF = w' * C_OFF_DIFF * w;

% Standard deviation lower bounds
sdGlxOFF = sqrt(varGlxOFF);
sdGlxJOINT = sqrt(varGlxJOINT);
sdGlxOFF_DIFF = sqrt(varGlxOFF_DIFF);

fprintf('\nGlx OFF SD: %.6f\n', sdGlxOFF);
fprintf('Glx JOINT Equal Time SD: %.6f\n', sdGlxJOINT);
fprintf('Glx OFF + DIFF SD: %.6f\n', sdGlxOFF_DIFF);

fprintf('\nJOINT / OFF: %.4f\n', ...
    sdGlxJOINT / sdGlxOFF);

fprintf('OFF + DIFF / OFF: %.4f\n', ...
    sdGlxOFF_DIFF / sdGlxOFF);