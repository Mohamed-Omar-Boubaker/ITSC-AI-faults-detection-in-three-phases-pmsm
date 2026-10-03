classdef ConsoleReport
    % All console printing of results (no computation)
    methods (Static)
        function hyperparameters(best, tl, cfg)
            fprintf('\n--- Selected Hyperparameters (RMSE-best among SV-ratio<=%.0f%% feasible trials) ---\n', cfg.TARGET_SV_RATIO*100);
            fprintf('  BoxConstraint = %.4f\n  KernelScale   = %.4f\n  Epsilon       = %.4f\n', ...
                best.BoxConstraint, best.KernelScale, best.Epsilon);
            ConsoleReport.nearBoundary('BoxConstraint', best.BoxConstraint, cfg.BOX_RANGE);
            ConsoleReport.nearBoundary('KernelScale',   best.KernelScale,   cfg.KS_RANGE);
            ConsoleReport.nearBoundary('Epsilon',       best.Epsilon,       cfg.EPS_RANGE);
            fprintf('\nTrials logged: %d | Feasible (SV<=%.0f%%) trials: %d\n', ...
                numel(tl), cfg.TARGET_SV_RATIO*100, sum([tl.sv_ratio] <= cfg.TARGET_SV_RATIO));
        end

        function nearBoundary(name, v, rg)
            span = rg(2) - rg(1);
            if (v - rg(1)) < 0.02*span || (rg(2) - v) < 0.02*span
                warning('%s = %.4f is within 2%% of its search range [%.2f, %.2f] -- consider widening this range.', ...
                    name, v, rg(1), rg(2));
            end
        end

        function regression(m, nTest)
            fprintf('\n==================================================\n REGRESSION PERFORMANCE METRICS (untouched test set)\n==================================================\n');
            fprintf('RMSE:      %.4f\nMAE:       %.4f\nMedAE:     %.4f\nMax Error: %.4f\n', m.rmse, m.mae, m.medae, m.maxerr);
            fprintf('MAPE:      %.2f%%  (%d/%d samples with y=0 excluded)\n', m.mape, sum(~m.mapeMask), nTest);
            fprintf('MAAPE:     %.2f%%  (all %d samples included)\n', m.maape, nTest);
            fprintf('R^2:       %.4f\n', m.r2);
        end

        function severity(s)
            fprintf('\n--- Severity-stratified RMSE (test set) ---\n');
            fprintf('Healthy RMSE: %.4f (n=%d) | Faulty RMSE: %.4f (n=%d)\n', s.healthy_rmse, s.healthy_n, s.faulty_rmse, s.faulty_n);
            fprintf('Low-sev RMSE: %.4f (n=%d) | High-sev RMSE: %.4f (n=%d)\n', s.low_rmse, s.low_n, s.high_rmse, s.high_n);
        end

        function binary(b, thr)
            fprintf('\n==================================================\n HEALTHY (0-%.2f%%) vs FAULTY (>%.2f%%) CLASSIFICATION (test set)\n==================================================\n', thr, thr);
            fprintf('Confusion Matrix:\n                    Predicted Healthy   Predicted Faulty\n');
            fprintf('Actual Healthy      TN = %-6d         FP = %-6d\n', b.TN, b.FP);
            fprintf('Actual Faulty       FN = %-6d         TP = %-6d\n\n', b.FN, b.TP);
            fprintf('Accuracy:              %.2f%%\nPrecision:             %.2f%%\nRecall (Sensitivity):  %.2f%%\n', b.accuracy, b.precision, b.recall);
            fprintf('Specificity:           %.2f%%\nF1 Score:              %.2f%%\nMatthews Corr (MCC):   %.2f%%\n', b.specificity, b.f1, b.mccPct);
        end

        function diagnostics(trRmse, rmse, trR2, r2, nSv, nTrain)
            fprintf('\n==================================================\n DIAGNOSTICS\n==================================================\n');
            fprintf('In-sample train RMSE: %.4f   |   Test RMSE: %.4f   |   Gap: %.4f\n', trRmse, rmse, rmse - trRmse);
            fprintf('In-sample train R^2:  %.4f   |   Test R^2:  %.4f   |   Gap: %.4f\n', trR2, r2, trR2 - r2);
            fprintf('Support vectors: %d / %d (%.1f%%)\n', nSv, nTrain, nSv/nTrain*100);
        end

        function perClass(titleStr, labels, m, C, width)
            fprintf('\n--- %s ---\n', titleStr);
            fmt = sprintf('  %%-%ds | Prec: %%5.1f%%%% | Rec: %%5.1f%%%% | Spec: %%5.1f%%%% | F1: %%5.1f%%%% | MCC: %%5.1f%%%% | N: %%d\n', width);
            for c = 1:numel(labels)
                fprintf(fmt, labels{c}, m.prec(c)*100, m.rec(c)*100, m.spec(c)*100, m.f1(c)*100, m.mcc(c)*100, sum(C(c, :)));
            end
        end

        function multiclassTest(testRes, sel)
            mf = testRes.metrics; K = sel.NClasses;
            fprintf('\n==================================================\n FINAL MULTICLASS TEST RESULTS (%d classes, frozen boundaries)\n==================================================\n', K);
            fprintf('Frozen boundaries: %s\n', mat2str(sel.Boundaries, 6));
            fprintf('Test Accuracy:          %.2f%%\nTest Balanced Accuracy: %.2f%%\n', testRes.accuracy, testRes.balanced_accuracy);
            fprintf('Test Macro F1:          %.2f%%\nTest Macro MCC:         %.2f%%\nWeighted F1:            %.2f%%\n', mf.macro_f1, mf.macro_mcc, mf.weight_f1);
            fprintf('\nConfusion Matrix (rows=true, cols=predicted):\n'); disp(testRes.C);
            ConsoleReport.perClass('Per-Class Test Metrics', sel.ClassNames, mf, testRes.C, 28);
            gap = sel.Chosen.accuracy - testRes.accuracy;
            fprintf('\nOOF Accuracy:  %.2f%%\nTest Accuracy: %.2f%%\nOOF-Test Gap:  %.2f%%\n', sel.Chosen.accuracy, testRes.accuracy, gap);
            if gap > 8, fprintf('  WARNING: possible boundary overfitting to OOF predictions.\n'); end
        end

        function targets(cfg, svRatio, bin, accMulti, K)
            fprintf('\n==================================================\n TARGET CHECK SUMMARY\n==================================================\n');
            fprintf('SV ratio target (<=%.0f%%):         %.1f%%   [%s]\n', cfg.TARGET_SV_RATIO*100, svRatio*100, ConsoleReport.passFail(svRatio <= cfg.TARGET_SV_RATIO));
            fprintf('Binary accuracy target (>=%.0f%%):   %.2f%%  [%s]\n', cfg.TARGET_ACCURACY*100, bin.accuracy, ConsoleReport.passFail(bin.accuracy/100 >= cfg.TARGET_ACCURACY));
            fprintf('%d-class accuracy target (>=%.0f%%): %.2f%%  [%s]  (min classes required: %d, max class size: %.0f%%, target class width: %.1f%%)\n', ...
                K, cfg.TARGET_ACCURACY*100, accMulti, ConsoleReport.passFail(accMulti/100 >= cfg.TARGET_ACCURACY), ...
                cfg.MIN_TOTAL_CLASSES, cfg.MAX_CLASS_FRAC*100, cfg.TARGET_CLASS_WIDTH);
        end

        function s = passFail(cond)
            if cond, s = 'PASS'; else, s = 'FAIL'; end
        end
    end
end
