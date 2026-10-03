classdef SvrSubsetEvaluator < SubsetEvaluator
    properties (SetAccess = immutable)
        C
        Epsilon
    end
    methods
        function obj = SvrSubsetEvaluator(nSplits, seed, C, epsilon)
            obj@SubsetEvaluator(nSplits, seed); obj.C = C; obj.Epsilon = epsilon;
        end
        function yhat = fitPredict(obj, Xtr, ytr, Xva)
            [mu, sd] = SvrFactory.scalerFit(Xtr);          % scaler fitted inside each fold
            A = (Xtr - mu) ./ sd; B = (Xva - mu) ./ sd;
            m = SvrFactory.fit(A, ytr, obj.C, obj.Epsilon, SvrFactory.gammaScale(A));
            yhat = predict(m, B);
        end
    end
end