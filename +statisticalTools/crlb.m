function [C, info] = crlb(F)
%CRLB Joint covariance lower bound from a real Fisher information matrix.
%   [C, info] = statisticalTools.crlb(F)
%
%   F : nParameters-by-nParameters real symmetric positive-semidefinite FIM.
%   C : inverse F if numerically full rank; NaN matrix otherwise.
%
%   info fields: rank, isFullRank, conditionNumber, reciprocalCondition,
%                tolerance, status.
%
%   This routine intentionally NEVER substitutes pinv(F) for a joint
%   covariance lower bound when F is singular. For an identifiable
%   combination of parameters even if F is singular, use
%   statisticalTools.contrastCRLB(F, L).

    if ~isnumeric(F) || ~ismatrix(F) || isempty(F) || ...
            size(F, 1) ~= size(F, 2) || ~isreal(F) || ...
            any(~isfinite(F(:)))
        error('statisticalTools:crlb:BadFIM', ...
            'F must be a finite, real, square matrix.');
    end
    asymmetry = norm(F - F.', 'fro') / max(norm(F, 'fro'), realmin);
    if asymmetry > 1e-10
        error('statisticalTools:crlb:AsymmetricFIM', ...
            'F must be symmetric.');
    end
    F = (F + F.') / 2;
    n = size(F, 1);
    eigValues = eig(F);
    tol = 10 * n * eps(max(abs(eigValues)));
    if any(eigValues < -tol)
        error('statisticalTools:crlb:IndefiniteFIM', ...
            'Fisher information must be positive semidefinite.');
    end
    r = nnz(eigValues > tol);
    info = struct('rank', r, 'isFullRank', r == n, ...
                  'conditionNumber', Inf, 'reciprocalCondition', rcond(F), ...
                  'tolerance', tol, 'status', 'rank-deficient');

    if r < n
        C = NaN(n);
        return
    end

    [R, flag] = chol(F);
    if flag ~= 0
        C = NaN(n);
        info.isFullRank = false;
        info.status = 'numerically-not-positive-definite';
        return
    end
    % F = R' * R, so inv(F) = inv(R) * inv(R').
    C = R \ (R' \ eye(n));
    C = (C + C.') / 2;
    info.conditionNumber = max(eigValues) / min(eigValues);
    info.status = 'ok';
end
