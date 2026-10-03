classdef ForestRandomSearch
    % Random search over forest settings scored by CV MSE, then refit on all data
    properties (SetAccess = immutable)
        NIter
        NSplits
        Seed
    end
    methods
        function obj = ForestRandomSearch(nIter, nSplits, seed)
            obj.NIter = nIter; obj.NSplits = nSplits; obj.Seed = seed;
        end
        function [model, params] = fit(obj, X, y, cols)
            rng(obj.Seed);
            depths = [Inf 8 12 20]; feats = {'sqrt', 'log2', 'all'};
            cvp = cvpartition(size(X, 1), 'KFold', obj.NSplits);
            best = Inf; params = [];
            for i = 1:obj.NIter
                q = struct('n_estimators', randi([100 349]), 'max_depth', depths(randi(4)), ...
                    'min_samples_split', randi([2 7]), 'min_samples_leaf', randi([1 5]), 'max_features', feats{randi(3)});
                e = zeros(cvp.NumTestSets, 1);
                for k = 1:cvp.NumTestSets
                    tr = training(cvp, k); va = test(cvp, k);
                    m = ForestFactory.fit(X(tr, cols), y(tr), q.n_estimators, q.max_depth, q.min_samples_split, q.min_samples_leaf, q.max_features);
                    e(k) = mean((y(va) - predict(m, X(va, cols))).^2);
                end
                mse = mean(e);
                fprintf('  search %d/%d  CV MSE = %.4f\n', i, obj.NIter, mse);
                if mse < best, best = mse; params = q; end
            end
            m = ForestFactory.fit(X(:, cols), y, params.n_estimators, params.max_depth, ...
                params.min_samples_split, params.min_samples_leaf, params.max_features);
            model = SelectedModel(m, cols, [], []);
        end
    end
end