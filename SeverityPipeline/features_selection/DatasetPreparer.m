classdef DatasetPreparer
    methods (Static)
        function [X, y] = clean(feat, percent)           % Inf -> NaN, drop NaN targets
            X = double(feat); y = double(percent(:));
            infMask = isinf(X);
            if any(infMask(:)), fprintf('Found %d Inf values -> converted to NaN\n', nnz(infMask)); X(infMask) = NaN; end
            nanCols = find(any(isnan(X), 1));
            if ~isempty(nanCols), fprintf('Columns with NaN values: %d (median-imputed using training data only)\n', numel(nanCols)); end
            bad = isnan(y);
            if any(bad), fprintf('Removing %d samples with NaN target\n', nnz(bad)); X = X(~bad, :); y = y(~bad); end
        end
        function [Xtr, Xte, ytr, yte] = splitAndImpute(X, y, seed)   % stratified 80/20 + train-median imputation
            s = rng; rng(seed);
            cvp = cvpartition(DatasetPreparer.quantileBins(y), 'HoldOut', 0.20);
            rng(s);
            tr = training(cvp); te = test(cvp);
            Xtr = X(tr, :); Xte = X(te, :); ytr = y(tr); yte = y(te);
            med = median(Xtr, 1, 'omitnan');
            for j = find(any(isnan(Xtr), 1) | any(isnan(Xte), 1))
                Xtr(isnan(Xtr(:, j)), j) = med(j); Xte(isnan(Xte(:, j)), j) = med(j);
            end
        end
        function bins = quantileBins(y)                  % like pandas.qcut (up to 10 bins, right-closed)
            nb = min(10, numel(unique(y))); bins = ones(size(y));
            if nb <= 1, return; end
            edges = unique(NumpyStats.percentile(y, linspace(0, 100, nb + 1)));
            for e = 2:numel(edges)-1, bins = bins + (y > edges(e)); end
        end
    end
end