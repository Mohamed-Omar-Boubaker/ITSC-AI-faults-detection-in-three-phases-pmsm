classdef SvrRandomSearch
    % Random search (log-uniform C, gamma, epsilon) scored by CV MSE, then refit on all data
    properties (SetAccess = immutable)
        NIter
        NSplits
        Seed
    end
    methods
        function obj = SvrRandomSearch(nIter, nSplits, seed)
            obj.NIter = nIter; obj.NSplits = nSplits; obj.Seed = seed;
        end
        function [model, params] = fit(obj, X, y, cols)
            [mu, sd] = SvrFactory.scalerFit(X(:, cols));
            Xs = (X(:, cols) - mu) ./ sd;
            rng(obj.Seed);
            u = @(lo, hi) exp(log(lo) + (log(hi) - log(lo)) * rand(obj.NIter, 1));
            C = u(1e-2, 1e3); gam = u(1e-4, 1e1); eps_ = u(1e-3, 1);
            cvp = cvpartition(size(Xs, 1), 'KFold', obj.NSplits);
            mse = zeros(obj.NIter, 1);
            for i = 1:obj.NIter
                e = zeros(cvp.NumTestSets, 1);
                for k = 1:cvp.NumTestSets
                    tr = training(cvp, k); va = test(cvp, k);
                    m = SvrFactory.fit(Xs(tr, :), y(tr), C(i), eps_(i), gam(i));
                    e(k) = mean((y(va) - predict(m, Xs(va, :))).^2);
                end
                mse(i) = mean(e);
                fprintf('  search %d/%d  CV MSE = %.4f\n', i, obj.NIter, mse(i));
            end
            [~, b] = min(mse);
            params = struct('C', C(b), 'gamma', gam(b), 'epsilon', eps_(b));
            model = SelectedModel(SvrFactory.fit(Xs, y, params.C, params.epsilon, params.gamma), cols, mu, sd);
        end
    end
end