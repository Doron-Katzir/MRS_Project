function Cout = covarianceTransform(C, M)
%COVARIANCETRANSFORM Propagate covariance through a real linear mapping.
%   Cout = statisticalTools.covarianceTransform(C, M)
%
%   Given real parameters theta with covariance C and derived quantities
%   z = M*theta, return Cov(z) = M*C*M'.
%
%   Example for theta = [Glu; Gln] and z = [Glu; Gln; Glx]:
%       M = [1 0; 0 1; 1 1];
%       Cz = statisticalTools.covarianceTransform(C, M);
%
%   For a SINGULAR Fisher matrix do not pass a pseudoinverse here;
%   use statisticalTools.contrastCRLB to identify estimable contrasts.

    if ~isnumeric(C) || ~ismatrix(C) || isempty(C) || ...
            ~isreal(C) || size(C,1) ~= size(C,2) || ...
            any(~isfinite(C(:))) || ...
            norm(C-C.', 'fro') > 1e-10*max(norm(C,'fro'),realmin)
        error('statisticalTools:covarianceTransform:BadCovariance', ...
            'C must be a finite, real, symmetric covariance matrix.');
    end
    if ~isnumeric(M) || ~ismatrix(M) || isempty(M) || ...
            ~isreal(M) || any(~isfinite(M(:))) || size(M,2) ~= size(C,1)
        error('statisticalTools:covarianceTransform:BadTransform', ...
            'M must have one column for every parameter in C.');
    end

    Cout = M * C * M.';
    Cout = (Cout + Cout.') / 2;
end
