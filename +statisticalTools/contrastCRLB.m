function [varianceBounds, info] = contrastCRLB(F, L)
%CONTRASTCRLB CRLB for real linear combinations of model parameters.
%   [varianceBounds, info] = statisticalTools.contrastCRLB(F, L)
%
%   F : real symmetric Fisher information (nParameters-by-nParameters).
%   L : nContrasts-by-nParameters. Each row defines a quantity L(i,:)*theta.
%
%   varianceBounds(i) is the variance lower bound for contrast i.
%   If F is singular, a contrast has a finite bound ONLY if its row
%   annihilates the null space of F; otherwise the bound is Inf.
%   A pseudoinverse is used ONLY for explicitly identifiable contrasts,
%   never as a replacement for the inverse joint FIM.
%
%   For two parameters [Glu; Gln], a useful contrast matrix is:
%       L = [1 0; 0 1; 1 1]; % Glu, Gln, Glx
%
%   info.identifiable is a logical vector and info.stdBounds gives SDs.

    if ~isnumeric(F) || ~ismatrix(F) || isempty(F) || ...
            ~isreal(F) || size(F, 1) ~= size(F, 2) || ...
            any(~isfinite(F(:))) || ...
            norm(F-F.', 'fro') > 1e-10*max(norm(F,'fro'),realmin)
        error('statisticalTools:contrastCRLB:BadFIM', ...
            'F must be finite, real, square and symmetric.');
    end
    n = size(F, 1);
    if ~isnumeric(L) || ~ismatrix(L) || ~isreal(L) || ...
            isempty(L) || size(L,2) ~= n || any(~isfinite(L(:)))
        error('statisticalTools:contrastCRLB:BadContrasts', ...
            'L must be a finite real matrix with one column per parameter.');
    end

    [V, D] = eig((F+F.') / 2);
    lambda = diag(D);
    tol = 10 * n * eps(max(abs(lambda)));
    if any(lambda < -tol)
        error('statisticalTools:contrastCRLB:IndefiniteFIM', ...
            'Fisher information must be positive semidefinite.');
    end
    good = lambda > tol;
    Vnull = V(:, ~good);
    Vgood = V(:, good);
    if any(good)
        Cidentifiable = (Vgood * diag(1 ./ lambda(good))) * Vgood.';
    else
        Cidentifiable = zeros(n);
    end

    nContrasts = size(L,1);
    varianceBounds = Inf(nContrasts, 1);
    identifiable = false(nContrasts, 1);
    for i = 1:nContrasts
        li = L(i,:);
        if norm(li * Vnull) <= 1e-8 * max(norm(li), realmin)
            identifiable(i) = true;
            varianceBounds(i) = max(0, li * Cidentifiable * li.');
        end
    end
    info = struct('identifiable', identifiable, ...
                  'stdBounds', sqrt(varianceBounds), ...
                  'informationRank', nnz(good), 'tolerance', tol);
end
