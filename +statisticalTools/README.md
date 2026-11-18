# `statisticalTools` — MATLAB package for MRS estimation

## Critical conventions

- `J` has **complex spectral samples in rows**, and **real-valued unknown parameters in columns**. For a known fixed linear model, the Jacobian with respect to concentrations is the basis matrix `B`.
- `noiseCov` is the **total complex noise variance** `E[abs(noise)^2]` (or covariance `E[noise*noise']`). For circular noise, the variance of the real channel is **half** that value.
- **Vertical concatenation** of independent ON/OFF observations uses `[J_ON; J_OFF]` with a diagonal variance vector or `blkdiag(C_ON,C_OFF)`. Concatenation **does not** amplify the noise SD of an individual observed point.
- **Subtraction** uses `J_DIFF = J_ON - J_OFF` and `C_DIFF = C_ON + C_OFF` for independent ON/OFF noise. If both have variance `v`, the DIFF variance is `2*v`.
- **Equal-time comparison** matters. If a total of `N` transients is available, OFF-only can average `N` transients, while an equal ON/OFF joint experiment has `N/2` transients per condition. Therefore each averaged ON/OFF *complex noise variance* is `2` times the OFF-only variance, and each noise SD is `sqrt(2)` times higher.
- `ReadBasis` produces **1024 samples from 512 acquired FID points** in the current simulation. The extra FFT frequency bins are *correlated*, not new independent observations. Use `fftNoiseCovariance` with indices matching the selected rows, **or** construct the native 512-point spectral grid and use its uncorrelated FFT noise. For an `fftshift`ed twofold zero-filled FFT, the native grid is `1:2:end`; confirm exported ordering before relying on this.
- `fftNoiseCovariance` assumes MATLAB's unnormalized FFT and iid **time-domain** complex noise. If the experimental spectrum is scaled, filtered, or apodized, the noise model must include the same processing. A physically calibrated **absolute** CRLB additionally requires consistent signal/basis amplitude units and noise units.
- These initial bounds treat the basis and linewidth/phase/frequency shift as **known**. Later, nuisance-parameter derivatives can be added as columns of `J` and marginalized by inverting the full Fisher matrix.
- The numerical CRLB assumes locally unconstrained real parameters, not boundary-constrained NNLS estimates.


## Minimal smoke test

```matlab
% Run from MRS_Project root after adding its root to the MATLAB path.
addpath(fullfile(pwd, 'MEGAPRESS_Stats_Package', 'tests')); % if extracted beside repo files
test_statisticalTools
```
