classdef NestedCrossValidator
    % Optional unbiased performance estimate (SvrConfig.RUN_NESTED_CV)
    properties (SetAccess = immutable)
        Cfg
        Weighter
    end
    methods
        function obj = NestedCrossValidator(cfg)
            obj.Cfg = cfg; obj.Weighter = SampleWeighter(cfg);
        end
        function run(obj, X, y)
            c = obj.Cfg;
            fprintf('\n================ NESTED CROSS-VALIDATION ================\n');
            rm = []; gap = []; svr = [];
            for rep = 1:c.NESTED_OUTER_REPEATS
                outer = cvpartition(numel(y), 'KFold', c.NESTED_OUTER_FOLDS);
                for k = 1:c.NESTED_OUTER_FOLDS
                    otr = outer.training(k); ova = outer.test(k);
                    sc = StandardScaler().fit(X(otr, :));
                    Xo = sc.transform(X(otr, :)); Xv = sc.transform(X(ova, :));
                    yo = y(otr); yv = y(ova); w = obj.Weighter.compute(yo);
                    inner = CrossValidation.repeatedKFold(numel(yo), c.NESTED_INNER_FOLDS, c.NESTED_INNER_REPEATS);
                    [pb, ~] = ConstrainedSvrSearch(c).run(Xo, yo, w, inner, c.NESTED_INNER_MAX_EVALS);
                    mdl = SvrTrainer.train(Xo, yo, pb, w);
                    mv = RegressionMetrics.compute(yv, SvrTrainer.estimate(mdl, Xv, true));
                    mt = RegressionMetrics.compute(yo, SvrTrainer.estimate(mdl, Xo, true));
                    svo = size(mdl.SupportVectors, 1) / numel(yo);
                    rm(end+1,1) = mv.rmse; gap(end+1,1) = mv.rmse - mt.rmse; svr(end+1,1) = svo; %#ok<AGROW>
                    fprintf('Outer fold  RMSE=%.4f  gap=%.4f  SVratio=%.1f%%  (BC=%.2f KS=%.2f Eps=%.3f)\n', ...
                        mv.rmse, gap(end), svo*100, pb.BoxConstraint, pb.KernelScale, pb.Epsilon);
                end
            end
            n = numel(rm); tcrit = tinv(0.975, n-1);
            fprintf('\n--- Nested CV Unbiased Estimate (n=%d outer folds) ---\n', n);
            fprintf('RMSE: %.4f ± %.4f (95%% CI)\n', mean(rm), tcrit*std(rm)/sqrt(n));
            fprintf('SV ratio: %.1f%% ± %.1f%%\n', mean(svr)*100, std(svr)*100);
        end
    end
end
