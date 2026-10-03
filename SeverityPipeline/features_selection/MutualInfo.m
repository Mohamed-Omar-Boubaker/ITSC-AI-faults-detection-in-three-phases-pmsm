classdef MutualInfo
    % k-NN (Kraskov) mutual information per column, same recipe as sklearn.mutual_info_regression
    methods (Static)
        function mi = regression(X, y, k, seed)
            s = rng; rng(seed);
            p = size(X, 2); mi = zeros(1, p);
            y = MutualInfo.prep(y(:));
            for j = 1:p, mi(j) = MutualInfo.pair(MutualInfo.prep(X(:, j)), y, k); end
            rng(s);
        end
        function v = prep(v)                      % scale by std, add tiny noise (breaks ties)
            sd = std(v, 1); if sd == 0, sd = 1; end
            v = v / sd;
            v = v + 1e-10 * max(1, mean(abs(v))) * randn(size(v));
        end
        function mi = pair(x, y, k)
            n = numel(x); Z = [x y];
            [~, D] = knnsearch(Z, Z, 'K', k + 1, 'Distance', 'chebychev');
            r  = D(:, end);                       % distance to the k-th neighbour
            nx = sum(abs(x - x') < r, 2) - 1;
            ny = sum(abs(y - y') < r, 2) - 1;
            mi = max(0, psi(n) + psi(k) - mean(psi(nx + 1)) - mean(psi(ny + 1)));
        end
    end
end