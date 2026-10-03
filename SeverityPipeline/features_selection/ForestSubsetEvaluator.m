classdef ForestSubsetEvaluator < SubsetEvaluator
    properties (SetAccess = immutable)
        Trees
        Depth
    end
    methods
        function obj = ForestSubsetEvaluator(nSplits, seed, trees, depth)
            obj@SubsetEvaluator(nSplits, seed); obj.Trees = trees; obj.Depth = depth;
        end
        function yhat = fitPredict(obj, Xtr, ytr, Xva)
            m = ForestFactory.fit(Xtr, ytr, obj.Trees, obj.Depth, 2, 1, 'all');
            yhat = predict(m, Xva);
        end
    end
end