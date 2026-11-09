%% Estimate the noise-limited fitting residual

% Choose a signal-free spectral region
% Verify visually that this is actually noise-only
noiseMask = ppm > 9 & ppm < 11;

noiseSamples = s(noiseMask);

% Estimate noise variance in each complex component
varReal = var(real(noiseSamples));
varImag = var(imag(noiseSamples));

% Expected complex noise power per frequency bin
noisePower = varReal + varImag;

% Number of complex points in the fitting window
K = numel(sFit);

% Expected residual norm from measurement noise
expectedNoiseNorm = sqrt(K * noisePower);

% Compare to measured fitting error
noiseBenchmark = 100 * expectedNoiseNorm / norm(sFit);

fprintf("Estimated noise benchmark: %.3f%%\n", ...
    noiseBenchmark);

fprintf("Actual residual: %.3f%%\n", ...
    100 * norm(rD));