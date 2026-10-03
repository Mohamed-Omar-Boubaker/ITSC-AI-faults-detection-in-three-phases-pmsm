classdef ConstrainedSvrSearch
    % Bayesian hyper-parameter search: CV-RMSE + penalty when the SV ratio exceeds the target
    properties (SetAccess = immutable)
        Cfg
    end
    methods
        function obj = ConstrainedSvrSearch(cfg), obj.Cfg = cfg; end

        function [best, trialLog] = run(obj, X, y, w, partitions, maxEvals)
            c = obj.Cfg; tl = TrialLog();
            vars = [optimizableVariable('BoxConstraint', c.BOX_RANGE, 'Transform', 'log'), ...
                    optimizableVariable('KernelScale',   c.KS_RANGE,  'Transform', 'log'), ...
                    optimizableVariable('Epsilon',       c.EPS_RANGE, 'Transform', 'log')];
            costFcn = @(p) obj.objective(p, X, y, w, partitions, tl);
            bayesopt(costFcn, vars, 'MaxObjectiveEvaluations', maxEvals, ...
                'IsObjectiveDeterministic', true, 'UseParallel', false, ...
                'AcquisitionFunctionName', 'expected-improvement-plus', ...
                'PlotFcn', [], 'Verbose', 1);
            trialLog = tl.Entries;
            best = obj.pickBestFeasible(trialLog);
        end

        function cost = objective(obj, params, X, y, w, partitions, tl)
            c = obj.Cfg; sqErr = []; totalSv = 0; totalN = 0;
            for r = 1:numel(partitions)
                cvp = partitions{r};
                for k = 1:cvp.NumTestSets
                    tr = cvp.training(k); va = cvp.test(k);
                    mdl = SvrTrainer.train(X(tr, :), y(tr), params, w(tr));
                    p = SvrTrainer.estimate(mdl, X(va, :), false);
                    sqErr = [sqErr; (y(va) - p).^2];                     %#ok<AGROW>
                    totalSv = totalSv + size(mdl.SupportVectors, 1);
                    totalN  = totalN + sum(tr);
                end
            end
            rmse = sqrt(mean(sqErr)); svRatio = totalSv / totalN;
            cost = rmse + c.SV_PENALTY_LAMBDA * max(0, svRatio - c.TARGET_SV_RATIO);
            tl.add(struct('BoxConstraint', params.BoxConstraint, 'KernelScale', params.KernelScale, ...
                'Epsilon', params.Epsilon, 'rmse', rmse, 'sv_ratio', svRatio, 'cost', cost));
            if c.VERBOSE_TRIALS
                if svRatio <= c.TARGET_SV_RATIO, feas = 'FEASIBLE'; else, feas = 'over-SV'; end
                fprintf('    [trial] RMSE=%.4f  SVratio=%.1f%% (%s)  BC=%.2f KS=%.2f Eps=%.3f  cost=%.4f\n', ...
                    rmse, svRatio*100, feas, params.BoxConstraint, params.KernelScale, params.Epsilon, cost);
            end
        end

        function p = pickBestFeasible(obj, trials)
            tgt = obj.Cfg.TARGET_SV_RATIO; ratios = [trials.sv_ratio];
            feasible = trials(ratios <= tgt);
            if isempty(feasible)
                [~, idx] = min(ratios);
                warning(['No trial met SV ratio <= %.0f%%. Falling back to the lowest-SV trial found (%.1f%%). ' ...
                    'Try: increasing MaxObjectiveEvaluations, widening BOX_RANGE/KS_RANGE, ' ...
                    'or reducing HEALTHY_WEIGHT_BOOST/HIGH_WEIGHT.'], tgt*100, trials(idx).sv_ratio*100);
                best = trials(idx);
            else
                [~, idx] = min([feasible.rmse]);
                best = feasible(idx);
            end
            p.BoxConstraint = best.BoxConstraint;
            p.KernelScale   = best.KernelScale;
            p.Epsilon       = best.Epsilon;
        end
    end
end
