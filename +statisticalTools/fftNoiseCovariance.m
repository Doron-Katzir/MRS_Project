function C = fftNoiseCovariance(nAcquired, nFFT, fftshiftIndices, timeNoiseVar)
%FFTNOISECOVARIANCE Noise covariance for FFT-shifted, zero-filled spectra.
%   C = statisticalTools.fftNoiseCovariance(nAcquired, nFFT, ...
%                                           fftshiftIndices, timeNoiseVar)
%
%   Assumes nAcquired independent proper complex time-domain noise samples
%   with E[abs(noise)^2] = timeNoiseVar. The remaining nFFT-nAcquired
%   samples are ZERO PADDING. The spectrum uses MATLAB's unnormalized FFT:
%       spectrum = fftshift(fft([fid; zeros(nFFT-nAcquired,1)]));
%
%   fftshiftIndices : indices into that full nFFT-point shifted spectrum.
%                     Use find(roiMask) if roiMask was built on that grid.
%   timeNoiseVar   : COMPLEX noise variance PER acquired time sample,
%                    AFTER any averaging. For N averages of iid noise,
%                    timeNoiseVar = singleAcquisitionVariance / N.
%
%   The output is Cov(spectralNoise) for the selected frequency bins.
%   Zero-filled spectral bins are CORRELATED. If more than nAcquired bins
%   are selected, the resulting covariance is necessarily singular;
%   statisticalTools.fisherInformation supports this if the corresponding
%   Jacobian lies in the same transform subspace.
%
%   Scaling and ordering must match the FFT used to create the spectra.
%   Any additional filtering, normalization, or apodization of measured
%   data must also be represented in its noise transformation.

    if ~isscalar(nAcquired) || ~isnumeric(nAcquired) || ...
            ~isreal(nAcquired) || ~isfinite(nAcquired) || ...
            nAcquired ~= floor(nAcquired) || nAcquired < 1 || ...
            ~isscalar(nFFT) || ~isnumeric(nFFT) || ~isreal(nFFT) || ...
            ~isfinite(nFFT) || nFFT ~= floor(nFFT) || nFFT < nAcquired
        error('statisticalTools:fftNoiseCovariance:BadLength', ...
            'Require integer nFFT >= nAcquired >= 1.');
    end
    if ~isnumeric(fftshiftIndices) || ~isreal(fftshiftIndices) || ...
            isempty(fftshiftIndices) || any(~isfinite(fftshiftIndices(:))) || ...
            any(fftshiftIndices(:) ~= floor(fftshiftIndices(:))) || ...
            any(fftshiftIndices(:) < 1) || any(fftshiftIndices(:) > nFFT) || ...
            numel(unique(fftshiftIndices(:))) ~= numel(fftshiftIndices)
        error('statisticalTools:fftNoiseCovariance:BadIndices', ...
            'Indices must be distinct integers between 1 and nFFT.');
    end
    if ~isscalar(timeNoiseVar) || ~isnumeric(timeNoiseVar) || ...
            ~isreal(timeNoiseVar) || ~isfinite(timeNoiseVar) || timeNoiseVar <= 0
        error('statisticalTools:fftNoiseCovariance:BadVariance', ...
            'timeNoiseVar must be a positive real scalar.');
    end

    % Each column of H is the spectrum resulting from a unit impulse at
    % one acquired time sample. Row selection matches the ROI selection.
    H = fftshift(fft(eye(nFFT, nAcquired), nFFT, 1), 1);
    H = H(fftshiftIndices(:), :);
    C = timeNoiseVar * (H * H');
    C = (C + C') / 2;
end
