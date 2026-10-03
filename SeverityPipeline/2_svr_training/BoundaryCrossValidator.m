classdef BoundaryCrossValidator
    % Honest CV score of the boundary search for a given number of classes
    properties (SetAccess = immutable)
        Cfg
        Optimizer
        Eval
    end
    methods
        function obj = BoundaryCrossValidator(cfg, opt, ev)
            obj.Cfg = cfg; obj.Optimizer = opt; obj.Eval = ev;
        end
        function s = score(obj, yTrain, oof, K, m)
            c = obj.Cfg; hb = obj.Eval.HealthyBound; N = numel(yTrain); scores = [];
            for rep = 1:c.BOUNDARY_CV_REPEATS
                cvp = cvpartition(N, 'KFold', c.BOUNDARY_CV_FOLDS);
                for k = 1:c.BOUNDARY_CV_FOLDS
                    fitI = cvp.training(k); valI = cvp.test(k);
                    ytf = yTrain(fitI); opf = oof(fitI); faultFit = ytf > hb;
                    if sum(faultFit) < 3 * max(m, 1), continue; end
                    minFit = max(5, round(c.MIN_CLASS_FRAC * sum(fitI)));
                    maxFit = round(c.MAX_CLASS_FRAC * sum(fitI));
                    [bf, found] = obj.Optimizer.optimize(ytf, opf, K, minFit, maxFit, max(opf(faultFit)));
                    if ~found, continue; end
                    minVal = max(3, round(c.MIN_CLASS_FRAC * sum(valI) * 0.5));
                    rv = obj.Eval.evaluate(yTrain(valI), oof(valI), bf.boundaries, K, minVal, inf);
                    if rv.feasible, scores(end+1) = rv.score; end        %#ok<AGROW>
                end
            end
            if isempty(scores), s = -Inf; else, s = mean(scores); end
        end
    end
end
