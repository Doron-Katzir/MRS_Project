classdef PostFitBasisCovariance
    %PostFitBasisCovariance Shared post-fit basis/covariance primitives.
    %
    % This class is the single interpretation of VDIIO.ReadLCMCoord
    % fitData.basisData used by both basis diagnostics. Reported sums are
    % excluded from the design matrix and propagated only after covariance
    % has been solved in the independent-component space.

    methods (Static)
        function [B, componentNames, info] = BuildIndependentBasis(fitData, quantRows)
            rawB = double(fitData.basisData);
            rawNames = string(fitData.basisMetName(:));

            if isempty(rawB) || size(rawB, 2) ~= numel(rawNames)
                error('fitData basisData/basisMetName dimensions are inconsistent.');
            end

            actualNames = strings(numel(rawNames), 1);
            concentration = nan(numel(rawNames), 1);
            keep = true(numel(rawNames), 1);
            reasons = strings(numel(rawNames), 1);
            quantNames = string(quantRows.name);

            for k = 1:numel(rawNames)
                [actualNames(k), qIdx] = ...
                    PostFitBasisCovariance.ResolveComponentName(rawNames(k), quantNames);

                if contains(actualNames(k), "+")
                    keep(k) = false;
                    reasons(k) = "reported sum is derived, not independent";
                    continue;
                end
                if isempty(qIdx)
                    keep(k) = false;
                    reasons(k) = "no matching fitted concentration";
                    continue;
                end

                concentration(k) = double(quantRows.sig(qIdx));
                if ~isfinite(concentration(k)) || concentration(k) == 0
                    keep(k) = false;
                    reasons(k) = "zero or nonfinite fitted concentration";
                elseif any(~isfinite(rawB(:, k)))
                    keep(k) = false;
                    reasons(k) = "nonfinite spectral samples";
                elseif norm(rawB(:, k)) == 0
                    keep(k) = false;
                    reasons(k) = "zero spectral column";
                end
            end

            if numel(unique(actualNames(keep))) ~= sum(keep)
                error('Duplicate active component names occur in fitted basis data.');
            end

            componentNames = actualNames(keep);
            componentConcentrations = concentration(keep);
            B = rawB(:, keep) ./ componentConcentrations.';

            excludedMask = ~keep;
            info = struct();
            info.actualComponentNames = actualNames;
            info.componentConcentrations = componentConcentrations;
            info.exclusionTable = table( ...
                rawNames(excludedMask), actualNames(excludedMask), reasons(excludedMask), ...
                'VariableNames', {'parsedName', 'resolvedName', 'reason'});
        end

        function [C, d] = SolveGeometry(designMatrix)
            designMatrix = double(designMatrix);
            H = real(designMatrix' * designMatrix);
            H = (H + H.') ./ 2;
            n = size(H, 1);

            d = struct();
            d.solverMethod = "unit-column-scaled thin QR";
            d.H = H;
            d.sizeH = size(H);
            d.rankH = rank(H);
            d.conditionNumber = cond(H);
            d.rcondH = rcond(H);
            d.rankA = rank(designMatrix);
            d.conditionNumberA = cond(designMatrix);

            eigenvalues = eig(H, 'vector');
            d.minimumEigenvalue = min(real(eigenvalues));
            d.maximumEigenvalue = max(real(eigenvalues));
            [~, cholStatus] = chol(H);
            d.positiveDefinite = cholStatus == 0;

            % Concentration normalization is scientific and has already been
            % applied by BuildIndependentBasis. This second scaling is purely
            % numerical: As = A*D^-1, D = diag(column L2 norms).
            columnScale = vecnorm(designMatrix, 2, 1);
            d.columnScale = columnScale;
            d.minimumColumnNorm = min(columnScale);
            d.maximumColumnNorm = max(columnScale);
            d.columnNormRatio = d.maximumColumnNorm / d.minimumColumnNorm;
            d.validColumnScale = all(isfinite(columnScale) & columnScale > 0);
            d.rankScaledA = NaN;
            d.conditionNumberScaledA = NaN;
            d.rcondScaledA = NaN;
            d.R = nan(n);
            d.absDiagonalR = nan(n, 1);
            d.smallestAbsDiagonalR = NaN;
            d.qrRankTolerance = NaN;
            d.rankFromQR = NaN;
            d.covarianceFinite = false;
            d.covariancePositiveDefinite = false;
            d.numericallyUsable = false;

            if ~d.validColumnScale || any(~isfinite(designMatrix), 'all')
                C = nan(n);
                warning('PostFitBasisCovariance:InvalidDesignScale', ...
                    ['Design matrix has a zero/nonfinite column norm or ', ...
                    'nonfinite entry. No pseudoinverse was used.']);
                return;
            end

            scaledDesign = designMatrix ./ columnScale;
            [~, R] = qr(scaledDesign, 0);
            absDiagonalR = abs(diag(R));
            qrTolerance = max(size(scaledDesign)) * eps(norm(R, inf));

            d.rankScaledA = rank(scaledDesign);
            d.conditionNumberScaledA = cond(scaledDesign);
            d.rcondScaledA = rcond(R);
            d.R = R;
            d.absDiagonalR = absDiagonalR;
            d.smallestAbsDiagonalR = min(absDiagonalR);
            d.qrRankTolerance = qrTolerance;
            d.rankFromQR = sum(absDiagonalR > qrTolerance);

            if d.rankFromQR ~= n || d.rankScaledA ~= n || ...
                    ~isfinite(d.conditionNumberScaledA)
                C = nan(n);
                warning('PostFitBasisCovariance:RankDeficientScaledDesign', ...
                    ['Scaled design matrix is numerically rank deficient: ', ...
                    'QR rank %d/%d, rank(As) %d/%d, cond(As) %.6g. ', ...
                    'No pseudoinverse was used.'], d.rankFromQR, n, ...
                    d.rankScaledA, n, d.conditionNumberScaledA);
                return;
            end

            % Cs = R^-1*R^-T, evaluated with a triangular solve. Because
            % As=A*D^-1 and beta_scaled=D*beta, transform covariance back as
            % C=D^-1*Cs*D^-1. No A'*A inverse is used for the solve.
            Rinverse = R \ eye(n);
            scaledCovariance = Rinverse * Rinverse.';
            inverseScale = 1 ./ columnScale;
            C = scaledCovariance .* (inverseScale.' * inverseScale);
            C = (C + C.') ./ 2;
            d.covarianceFinite = all(isfinite(C), 'all');
            [~, covarianceCholStatus] = chol(C);
            d.covariancePositiveDefinite = covarianceCholStatus == 0;
            d.numericallyUsable = d.covarianceFinite;

            if ~d.numericallyUsable
                C = nan(n);
                warning('PostFitBasisCovariance:NonfiniteCovariance', ...
                    'Scaled-QR covariance contains nonfinite values.');
            end
        end

        function [M, propagationTable] = BuildReportedTransformation( ...
                componentNames, reportedNames, candidateNames)
            nReported = numel(reportedNames);
            nComponents = numel(componentNames);
            M = zeros(nReported, nComponents);
            termsText = strings(nReported, 1);
            activeTermsText = strings(nReported, 1);
            inactiveTermsText = strings(nReported, 1);

            for r = 1:nReported
                terms = string(split(reportedNames(r), "+"));
                terms = terms(:);
                termsText(r) = strjoin(terms, " + ");
                activeTerms = strings(0, 1);
                inactiveTerms = strings(0, 1);

                for t = 1:numel(terms)
                    idx = find(componentNames == terms(t), 1, 'first');
                    if ~isempty(idx)
                        M(r, idx) = M(r, idx) + 1;
                        activeTerms(end+1, 1) = terms(t); %#ok<AGROW>
                    else
                        if ~ismember(terms(t), candidateNames)
                            error('Reported term %s is not a known independent candidate.', terms(t));
                        end
                        inactiveTerms(end+1, 1) = terms(t); %#ok<AGROW>
                    end
                end
                activeTermsText(r) = strjoin(activeTerms, " + ");
                inactiveTermsText(r) = strjoin(inactiveTerms, " + ");
            end

            propagationTable = table( ...
                reportedNames, termsText, activeTermsText, inactiveTermsText, ...
                sum(M ~= 0, 2), ...
                'VariableNames', {'reportedName', 'definedTerms', 'activeTerms', ...
                'inactiveFixedTerms', 'nActiveTerms'});
        end

        function baseline = GetBaselineColumn(fitData, spectralPointCount)
            baseline = double(fitData.baseline(:));
            if numel(baseline) ~= spectralPointCount
                error('Fitted baseline length (%d) does not match basis rows (%d).', ...
                    numel(baseline), spectralPointCount);
            end
            if any(~isfinite(baseline))
                error('Fitted baseline contains nonfinite spectral samples.');
            end
            if norm(baseline) == 0
                error('Fitted baseline is a zero spectral column.');
            end
        end

        function R = CovarianceToCorrelation(C)
            C = double(C);
            C = (C + C.') ./ 2;
            variance = diag(C);
            denom = sqrt(variance * variance.');
            R = C ./ denom;
            invalid = ~isfinite(variance) | variance <= 0;
            R(invalid, :) = NaN;
            R(:, invalid) = NaN;
            valid = ~invalid;
            diagonalIndices = find(valid);
            R(sub2ind(size(R), diagonalIndices, diagonalIndices)) = 1;
            R = (R + R.') ./ 2;

            finiteMask = isfinite(R);
            tolerance = 1e-10;
            if any(abs(R(finiteMask)) > 1 + tolerance)
                warning('A covariance-derived correlation exceeded [-1,1] beyond tolerance.');
            end
            R(finiteMask) = min(1, max(-1, R(finiteMask)));
        end

        function [fitData, quantRows] = GetLoadedPartData(covPatient, coordFile)
            coordFiles = string(covPatient.coordFiles(:));
            targetFull = lower(replace(string(coordFile), '/', '\'));
            normalized = lower(replace(coordFiles, '/', '\'));
            idx = find(normalized == targetFull, 1, 'first');

            if isempty(idx)
                targetName = PostFitBasisCovariance.GetFileName(coordFile);
                names = arrayfun(@(x) PostFitBasisCovariance.GetFileName(x), coordFiles);
                idx = find(strcmpi(names, targetName), 1, 'first');
            end
            if isempty(idx)
                error('Loaded fitData does not contain %s.', coordFile);
            end

            fitData = covPatient.fitData(idx);
            T = covPatient.coordTable;
            targetName = PostFitBasisCovariance.GetFileName(coordFile);
            rowMask = strcmpi(string(T.filename), targetName);
            quantRows = T(rowMask, :);
            if isempty(quantRows)
                error('Loaded quantification table has no rows for %s.', targetName);
            end
        end

        function names = DefaultReportedNames()
            names = [ ...
                "NAA"; "NAAG"; "Cr"; "PCr"; "GPC"; "PCh"; "Glu"; "Gln"; ...
                "GABA"; "GSH"; "Tau"; "Asc"; "Glc"; "Ace"; "mI"; "sI"; ...
                "Asp"; "Lac"; "GPC+PCh"; "NAA+NAAG"; "Cr+PCr"; "Glu+Gln"];
        end

        function names = DefaultIndependentCandidateNames()
            names = [ ...
                "NAA"; "NAAG"; "Cr"; "PCr"; "GPC"; "PCh"; "Glu"; "Gln"; ...
                "GABA"; "GSH"; "Tau"; "Asc"; "Glc"; "Ace"; "mI"; "sI"; ...
                "Asp"; "Lac"; "Lip13a"; "Lip13b"; "Lip09"; "MM09"; ...
                "Lip20"; "MM20"; "MM12"; "MM14"; "MM17"; "-CrCH2"];
        end
    end

    methods (Static, Access = private)
        function [resolved, idx] = ResolveComponentName(parsedName, quantNames)
            parsedName = string(parsedName);
            idx = find(strcmpi(quantNames, parsedName), 1, 'first');
            resolved = parsedName;

            % VDIIO's basis-header expression does not retain the leading
            % minus in LCModel's -CrCH2 label.
            if isempty(idx) && strcmpi(parsedName, "CrCH2")
                idx = find(strcmpi(quantNames, "-CrCH2"), 1, 'first');
                if ~isempty(idx)
                    resolved = "-CrCH2";
                end
            end
        end

        function name = GetFileName(pathValue)
            [~, base, ext] = fileparts(char(pathValue));
            name = string([base, ext]);
        end
    end
end
