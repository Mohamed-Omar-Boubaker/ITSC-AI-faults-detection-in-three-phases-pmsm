classdef ModelDiagnostics
    % Sanity / robustness checks. New check = new method + one line in runAll().
    properties (SetAccess = immutable)
        Cfg
        Model
        Params
        Weighter
        Xtr
        ytr
        Xte
        yte
        nTotal
        nUnique
    end
    methods
        function obj = ModelDiagnostics(cfg, model, params, weighter, Xtr, ytr, Xte, yte, yAll)
            obj.Cfg = cfg; obj.Model = model; obj.Params = params; obj.Weighter = weighter;
            obj.Xtr = Xtr; obj.ytr = ytr; obj.Xte = Xte; obj.yte = yte;
            obj.nTotal = numel(yAll); obj.nUnique = numel(unique(yAll));
        end

        function runAll(obj, m, yPred)
            obj.leakageChecks(m);
            obj.svComposition();
            obj.noiseRobustness(yPred);
            obj.learningCurve();
            obj.calibration(yPred);
        end

        function leakageChecks(obj, m)
            yS = obj.ytr(randperm(numel(obj.ytr)));
            mdl = SvrTrainer.train(obj.Xtr, yS, obj.Params, obj.Weighter.compute(yS));
            p = SvrTrainer.estimate(mdl, obj.Xte, true);
            rmseS = sqrt(mean((obj.yte - p).^2));
            r2S = 1 - sum((obj.yte - p).^2) / sum((obj.yte - mean(obj.yte)).^2);
            fprintf('\nShuffled-label sanity check: RMSE=%.4f (real=%.4f), R2=%.4f (real=%.4f)\n', rmseS, m.rmse, r2S, m.r2);
            if r2S > 0.3, fprintf('  WARNING: shuffled-label R^2 is suspiciously high -> possible leakage.\n'); end
            md = min(pdist2(obj.Xte, obj.Xtr, 'euclidean'), [], 2);
            fprintf('Test-to-train nearest neighbor distances: min=%.4f median=%.4f\n', min(md), median(md));
            fprintf('%d / %d test samples have a near-duplicate in training (dist<0.3)\n', sum(md < 0.3), numel(md));
            fprintf('Total samples: %d | Unique percent values: %d\n', obj.nTotal, obj.nUnique);
        end

        function svComposition(obj)
            fprintf('\n==================================================\n SUPPORT VECTOR COMPOSITION ANALYSIS\n==================================================\n');
            try
                mdl = obj.Model; svMask = mdl.IsSupportVector; alphaAbs = abs(mdl.Alpha);
                boxSv = mdl.BoxConstraints(svMask);
                if numel(boxSv) == numel(alphaAbs)
                    fb = mean(alphaAbs >= (boxSv - 1e-6)) * 100;
                    fprintf('Support vectors: %d total | Bounded (alpha=C, tube-violating): %.1f%% | Free (0<alpha<C, on-margin): %.1f%%\n', numel(alphaAbs), fb, 100 - fb);
                    if fb > 70, fprintf('  NOTE: most SVs are bounded -> model is straining to fit; consider larger Epsilon or lower sample weights.\n'); end
                else
                    fprintf('Could not align BoxConstraints with Alpha (sizes: %d vs %d) -- skipping bounded/free split.\n', numel(boxSv), numel(alphaAbs));
                end
                y = obj.ytr; thr = obj.Cfg.SEVERITY_THRESHOLD;
                coarse = double(y == 0) + 2*double(y > 0 & y <= thr) + 3*double(y > thr);
                lbl = {'Healthy (y=0)', 'Low fault (0<y<=thresh)', 'High fault (y>thresh)'};
                fprintf('\nPer coarse-class SV membership (train set):\n');
                for l = 1:3
                    mk = coarse == l; n = sum(mk); ns = sum(mk & svMask);
                    fprintf('  %-26s | N=%4d | SVs=%4d (%.1f%% of class became SVs)\n', lbl{l}, n, ns, 100*ns/max(n, 1));
                end
            catch ME
                fprintf('SV composition analysis skipped (property unavailable in this MATLAB version): %s\n', ME.message);
            end
        end

        function noiseRobustness(obj, yPred)
            fprintf('\n==================================================\n PREDICTION ROBUSTNESS TO FEATURE NOISE (test set)\n==================================================\n');
            thr = obj.Cfg.HEALTHY_THRESHOLD;
            for nl = [0.02, 0.05, 0.10]
                yn = SvrTrainer.estimate(obj.Model, obj.Xte + nl * randn(size(obj.Xte)), true);
                flips = sum((yn > thr) ~= (yPred > thr));
                fprintf('Noise sigma=%.2f: mean prediction shift RMS=%.4f | healthy/faulty decision flips=%d/%d (%.1f%%)\n', ...
                    nl, sqrt(mean((yn - yPred).^2)), flips, numel(yPred), 100*flips/numel(yPred));
            end
        end

        function learningCurve(obj)
            fprintf('\n==================================================\n LEARNING CURVE (fixed test set, growing training size)\n==================================================\n');
            fr = [0.2, 0.4, 0.6, 0.8, 1.0]; n = zeros(size(fr)); rm = n; sv = n;
            fprintf('%-10s %-12s %-12s %-10s\n', 'TrainFrac', 'TrainN', 'TestRMSE', 'SVratio%');
            for i = 1:numel(fr)
                nSub = round(fr(i) * numel(obj.ytr)); idx = randperm(numel(obj.ytr), nSub);
                ys = obj.ytr(idx);
                mdl = SvrTrainer.train(obj.Xtr(idx, :), ys, obj.Params, obj.Weighter.compute(ys));
                p = SvrTrainer.estimate(mdl, obj.Xte, true);
                n(i) = nSub; rm(i) = sqrt(mean((obj.yte - p).^2)); sv(i) = size(mdl.SupportVectors, 1) / nSub * 100;
                fprintf('%-10.0f%% %-12d %-12.4f %-10.1f\n', fr(i)*100, nSub, rm(i), sv(i));
            end
            ResultPlotter.learningCurve(n, rm, sv);
        end

        function calibration(obj, yPred)
            fprintf('\n==================================================\n CALIBRATION CHECK: predicted vs actual severity, by bin (test set)\n==================================================\n');
            nb = 5; [~, edges] = histcounts(yPred, nb); bi = discretize(yPred, edges);
            fprintf('%-8s %-12s %-14s %-10s\n', 'Bin', 'MeanPred', 'MeanActual', 'N');
            for b = 1:nb
                mk = bi == b;
                if any(mk), fprintf('%-8d %-12.3f %-14.3f %-10d\n', b, mean(yPred(mk)), mean(obj.yte(mk)), sum(mk)); end
            end
        end
    end
end
