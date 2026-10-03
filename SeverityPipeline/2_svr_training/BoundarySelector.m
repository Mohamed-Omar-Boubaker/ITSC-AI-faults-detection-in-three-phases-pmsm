classdef BoundarySelector
    % Picks the number of classes (complexity-penalised CV score) or falls back to binary
    properties (SetAccess = immutable)
        Cfg
        Optimizer
        CvScorer
    end
    methods
        function obj = BoundarySelector(cfg, opt, cvScorer)
            obj.Cfg = cfg; obj.Optimizer = opt; obj.CvScorer = cvScorer;
        end

        function sel = select(obj, yTrain, oof, bin)
            c = obj.Cfg; hb = c.HEALTHY_THRESHOLD; N = numel(yTrain);
            maxFaultPred = max(oof(yTrain > hb));
            minSize = max(15, round(c.MIN_CLASS_FRAC * N)); maxSize = round(c.MAX_CLASS_FRAC * N);
            fprintf('Class size constraint: %d <= class N <= %d (out of %d train samples)\n', minSize, maxSize, N);
            fprintf('Target class width: %.1f%% (soft, bonus weight=%.2f)\n', c.TARGET_CLASS_WIDTH, c.WIDTH_BONUS_WEIGHT);

            allK = struct('K', {}, 'boundaries', {}, 'res', {}, 'cv_score', {});
            for K = c.MIN_TOTAL_CLASSES:c.MAX_TOTAL_CLASSES
                [res, found] = obj.Optimizer.optimize(yTrain, oof, K, minSize, maxSize, maxFaultPred);
                if found
                    cvs = obj.CvScorer.score(yTrain, oof, K, K - 2);
                    fprintf('  %d classes: FullFitScore=%.4f  HonestCVScore=%.4f  boundaries=%s\n', ...
                        K, res.score, cvs, mat2str(res.boundaries, 6));
                    allK(end+1) = struct('K', K, 'boundaries', res.boundaries, 'res', res, 'cv_score', cvs); %#ok<AGROW>
                else
                    fprintf('  %d classes: NO VALID SOLUTION\n', K);
                end
            end

            if isempty(allK)                                   % ---- binary fallback ----
                fprintf('\nNo K in [%d,%d] is feasible -- falling back to binary Healthy/Faulty.\n', c.MIN_TOTAL_CLASSES, c.MAX_TOTAL_CLASSES);
                sel.NClasses = 2; sel.Edges = [-inf, hb, inf]; sel.Boundaries = [];
                sel.ClassNames = {sprintf('Healthy (0-%.2f%%)', hb), sprintf('Faulty (>%.2f%%)', hb)};
                sel.Chosen = struct('n_total_classes', 2, 'edges', sel.Edges, 'accuracy', bin.accuracy, ...
                    'macro_f1', bin.f1, 'macro_mcc', bin.mccPct);
                sel.Boot = struct('mean_acc', NaN, 'std_acc', NaN, 'ci_lower', NaN, 'ci_upper', NaN);
                return;
            end

            adj = zeros(1, numel(allK));
            fprintf('\n--- Complexity-penalized comparison (lambda=%.4f per extra class) ---\n', c.COMPLEXITY_LAMBDA);
            for i = 1:numel(allK)
                pen = c.COMPLEXITY_LAMBDA * (allK(i).K - c.MIN_TOTAL_CLASSES);
                adj(i) = allK(i).cv_score - pen;
                fprintf('  %d classes: CVScore=%.4f - Penalty=%.4f = %.4f\n', allK(i).K, allK(i).cv_score, pen, adj(i));
            end
            [~, bi] = max(adj); g = allK(bi);
            fprintf('\n>>> SELECTED: %d total classes | OOF Accuracy=%.2f%%\n', g.K, g.res.accuracy);

            sel.Boundaries = g.boundaries; sel.Edges = g.res.edges; sel.NClasses = g.K;
            sel.ClassNames = BoundarySelector.classNames(sel.Edges, hb);
            fprintf('\nFrozen class boundaries:\n'); fprintf('  %s\n', sel.ClassNames{:});
            sel.Chosen = struct('n_total_classes', g.K, 'edges', sel.Edges, 'accuracy', g.res.accuracy, ...
                'macro_f1', g.res.macro_f1, 'macro_mcc', g.res.macro_mcc);
            sel.Boot = BootstrapStability.run(yTrain, oof, sel.Edges, c.N_BOOT);
            fprintf('\nBootstrap (OOF, %d resamples): Accuracy = %.2f%% +/- %.2f%% (95%% CI [%.2f%%, %.2f%%])\n', ...
                c.N_BOOT, sel.Boot.mean_acc, sel.Boot.std_acc, sel.Boot.ci_lower, sel.Boot.ci_upper);
            if (sel.Boot.ci_upper - sel.Boot.ci_lower) > 10
                fprintf('  NOTE: wide bootstrap CI -- boundary placement is sensitive to which samples land in training.\n');
            end
        end
    end
    methods (Static)
        function names = classNames(edges, hb)
            K = numel(edges) - 1; names = cell(1, K);
            names{1} = sprintf('Healthy (0-%.2f%%)', hb);
            for c = 2:K
                lo = edges(c); hi = edges(c+1);
                if isinf(hi), names{c} = sprintf('Class%d (>%.2f%%)', c, lo);
                else,         names{c} = sprintf('Class%d (%.2f-%.2f%%)', c, lo, hi);
                end
            end
        end
    end
end
