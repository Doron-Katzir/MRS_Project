function [F, info] = fisherInformation(J, noiseCov)
%FISHERINFORMATION Fisher information for real parameters, complex data.
%   F = statisticalTools.fisherInformation(J, noiseCov)
%   [F, info] = statisticalTools.fisherInformation(J, noiseCov)
%
%   Model: y = mu(theta) + noise, with theta REAL and noise a proper
%   (circular) complex Gaussian vector whose known covariance does not
%   depend on theta.
%
%   J        : nObservations-by-nParameters complex Jacobian d(mu)/d(theta).
%   noiseCov : either
%                - positive scalar (equal complex variance per observation),
%                - positive nObservations-vector (independent unequal noise),
%                - nObservations-by-nObservations Hermitian covariance.
%
%   Convention: noiseCov = E[noise * noise'], NOT the variance of just
%   the real channel. For circular noise, each real/imaginary channel has
%   half of this variance.
%
%   F = 2 * real(J' * inv(noiseCov) * J), evaluated by whitening rather
%   than explicitly inverting the covariance matrix.
%
%   A singular (positive-semidefinite) full covariance is allowed only
%   when EVERY column of J lies in the covariance's column space. This
%   is useful for zero-filled FFT spectra. An inconsistent J causes an
%   error instead of pretending a zero-noise direction is informative.
%
%   Assumptions: proper complex Gaussian noise; real unknown parameters;
%   known, parameter-independent covariance. Do not use this formula
%   unchanged for improper noise or covariance-dependent parameters.

    if ~isnumeric(J) || ~ismatrix(J) || isempty(J) || ...
            any(~isfinite(J(:)))
        error('statisticalTools:fisherInformation:BadJacobian', ...
            'J must be a finite, nonempty numeric 2-D matrix.');
    end

    nObs = size(J, 1);
    info = struct('noiseModel', '', 'noiseRank', nObs, ...
                  'numObservations', nObs, ...
                  'numParameters', size(J, 2));

    if ~isnumeric(noiseCov) || isempty(noiseCov) || ...
            any(~isfinite(noiseCov(:)))
        error('statisticalTools:fisherInformation:BadNoise', ...
            'Noise covariance must be finite and numeric.');
    end

    if isscalar(noiseCov)
        if ~isreal(noiseCov) || noiseCov <= 0
            error('statisticalTools:fisherInformation:BadNoise', ...
                'Complex noise variance must be positive and real.');
        end
        Jw = J / sqrt(noiseCov);
        info.noiseModel = 'iid';

    elseif isvector(noiseCov)
        v = noiseCov(:);
        if numel(v) ~= nObs || ~isreal(v) || any(v <= 0)
            error('statisticalTools:fisherInformation:BadNoiseVector', ...
                'Independent-noise variance vector must have one positive real entry per observation.');
        end
        Jw = bsxfun(@rdivide, J, sqrt(v));
        info.noiseModel = 'independent-unequal';

    else
        C = noiseCov;
        if ~ismatrix(C) || ~isequal(size(C), [nObs, nObs])
            error('statisticalTools:fisherInformation:BadNoiseSize', ...
                'Full noise covariance must be nObservations-by-nObservations.');
        end
        hermError = norm(C - C', 'fro') / max(norm(C, 'fro'), realmin);
        if hermError > 1e-10
            error('statisticalTools:fisherInformation:NonHermitian', ...
                'Full noise covariance must be Hermitian.');
        end
        C = (C + C') / 2;
        info.noiseModel = 'full-covariance';

        [L, cholFlag] = chol(C, 'lower');
        % A theoretically singular covariance can occasionally pass chol
        % because roundoff gives its zero eigenvalues tiny positive pivots.
        % Detect near singularity before trusting a Cholesky whitening.
        if cholFlag == 0 && rcond(C) > 1e-12
            Jw = L \ J;
        else
            % The covariance can be singular after zero-filling an FFT.
            % Whiten only the stochastic subspace, and verify the model
            % has no derivative outside that subspace.
            [U, D] = eig(C);
            lambda = real(diag(D));
            tol = 10 * nObs * eps(max(abs(lambda)));
            if any(lambda < -tol)
                error('statisticalTools:fisherInformation:NotPSD', ...
                    'Noise covariance must be positive semidefinite.');
            end
            good = lambda > tol;
            info.noiseRank = nnz(good);
            if ~any(good)
                error('statisticalTools:fisherInformation:NoNoiseSubspace', ...
                    'Noise covariance contains no positive-variance directions.');
            end
            nullProjection = U(:, ~good)' * J;
            relativeMismatch = norm(nullProjection, 'fro') / ...
                max(norm(J, 'fro'), realmin);
            if relativeMismatch > 1e-7
                error('statisticalTools:fisherInformation:SignalOutsideNoiseSubspace', ...
                    ['J has components outside the noise covariance range. ', ...
                     'Check FFT ordering, scaling, or the noise model.']);
            end
            Jw = bsxfun(@rdivide, U(:, good)' * J, sqrt(lambda(good)));
        end
    end

    F = 2 * real(Jw' * Jw);
    F = (F + F.') / 2;
end
