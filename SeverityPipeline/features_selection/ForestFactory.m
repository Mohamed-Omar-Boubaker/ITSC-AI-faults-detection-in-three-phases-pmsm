classdef ForestFactory
    % sklearn-like RandomForestRegressor. MATLAB has no max_depth: depth d -> at most 2^d-1 splits.
    methods (Static)
        function m = fit(X, y, nTrees, maxDepth, minSplit, minLeaf, maxFeat)
            p = size(X, 2); n = size(X, 1);
            switch maxFeat
                case 'sqrt', nv = max(1, floor(sqrt(p)));
                case 'log2', nv = max(1, floor(log2(p)));
                otherwise,   nv = p;
            end
            if isinf(maxDepth), ms = n - 1; else, ms = min(2^maxDepth - 1, n - 1); end
            t = templateTree('MaxNumSplits', ms, 'MinParentSize', minSplit, ...
                'MinLeafSize', minLeaf, 'NumVariablesToSample', nv);
            m = fitrensemble(X, y, 'Method', 'Bag', 'NumLearningCycles', nTrees, 'Learners', t);
        end
    end
end