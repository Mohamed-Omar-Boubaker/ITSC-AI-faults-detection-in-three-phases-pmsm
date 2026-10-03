classdef BoundaryEvaluator
    % Scores one set of class boundaries (the healthy bound is fixed, not searched)
    properties (SetAccess = immutable)
        HealthyBound
        ScoreWeights
        UseBalanced
    end
    methods
        function obj = BoundaryEvaluator(hb, w, useBal)
            obj.HealthyBound = hb; obj.ScoreWeights = w; obj.UseBalanced = useBal;
        end
        function res = evaluate(obj, yTrue, yHat, mid, K, minSize, maxSize)
            edges = [-inf, obj.HealthyBound, mid, inf];
            tc = discretize(yTrue, edges); pc = discretize(yHat, edges);
            res = struct('feasible', false, 'score', -Inf, 'accuracy', NaN, 'balanced_accuracy', NaN, ...
                'macro_f1', NaN, 'macro_mcc', NaN, 'weight_f1', NaN, 'counts', [], 'C', [], ...
                'metrics', [], 'edges', edges);
            if any(isnan(tc)) || any(isnan(pc)), return; end
            counts = accumarray(tc, 1, [K, 1]); res.counts = counts;
            fc = counts(2:end);
            res.feasible = ~any(fc < minSize) && ~any(fc > maxSize);
            if ~res.feasible, return; end
            C = confusionmat(tc, pc, 'Order', 1:K);
            m = ConfusionMetrics.multiClass(C);
            acc = mean(tc == pc); bal = ConfusionMetrics.balancedAccuracy(C);
            if obj.UseBalanced, accTerm = bal; else, accTerm = acc; end
            w = obj.ScoreWeights;
            res.accuracy = acc * 100; res.balanced_accuracy = bal * 100;
            res.macro_f1 = m.macro_f1; res.macro_mcc = m.macro_mcc; res.weight_f1 = m.weight_f1;
            res.C = C; res.metrics = m;
            res.score = w(1)*accTerm + w(2)*(m.macro_f1/100) + w(3)*(m.macro_mcc/100);
        end
    end
end
