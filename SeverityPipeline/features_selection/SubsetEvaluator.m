classdef (Abstract) SubsetEvaluator
    % Cross-validated RMSE of a column subset. Same folds for every call (fixed seed).
    properties (SetAccess = immutable)
        NSplits
        Seed
    end
    methods (Abstract)
        yhat = fitPredict(obj, Xtr, ytr, Xva)
    end
    methods
        function obj = SubsetEvaluator(nSplits, seed), obj.NSplits = nSplits; obj.Seed = seed; end
        function r = cvRmse(obj, X, y, cols)
            s = rng; rng(obj.Seed); cvp = cvpartition(size(X, 1), 'KFold', obj.NSplits); rng(s);
            yp = zeros(size(y));
            for k = 1:cvp.NumTestSets
                tr = training(cvp, k); va = test(cvp, k);
                yp(va) = obj.fitPredict(X(tr, cols), y(tr), X(va, cols));
            end
            r = sqrt(mean((y - yp).^2));
        end
    end
end