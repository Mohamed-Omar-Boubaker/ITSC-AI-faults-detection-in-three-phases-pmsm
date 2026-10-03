classdef ResultPlotter
    % One static method per figure. R = results struct built by SeverityExperiment.
    methods (Static)
        function all(R)
            set(0, 'DefaultFigureColor', [1 1 1]);
            ResultPlotter.predictedVsActual(R);
            ResultPlotter.residuals(R);
            ResultPlotter.orderedTracking(R);
            ResultPlotter.cvFolds(R);
            ResultPlotter.generalizationGap(R);
            ResultPlotter.permutationImportance(R);
            ResultPlotter.searchHistory(R);
            ResultPlotter.multiclassDiagnostics(R);
            ResultPlotter.perClassMetrics(R);
            ResultPlotter.rocCurves(R);
        end

        function setLabels(xl, yl, tt)
            xlabel(xl, 'FontSize', 11, 'FontWeight', 'bold');
            ylabel(yl, 'FontSize', 11, 'FontWeight', 'bold');
            title(tt, 'FontSize', 12);
        end

        function learningCurve(n, rm, sv)
            figure('Name', 'Learning Curve - Test RMSE');
            plot(n, rm, '-o', 'LineWidth', 1.5, 'MarkerFaceColor', 'auto');
            xlabel('Training Set Size (samples)'); ylabel('Test RMSE');
            title('Learning Curve: Test RMSE vs. Training Set Size'); grid on;
            figure('Name', 'Learning Curve - SV Ratio');
            plot(n, sv, '-s', 'LineWidth', 1.5, 'Color', [0.85 0.33 0.10], 'MarkerFaceColor', 'auto');
            xlabel('Training Set Size (samples)'); ylabel('Support Vector Ratio (%)');
            title('Learning Curve: SV Ratio vs. Training Set Size'); grid on;
        end

        function predictedVsActual(R)
            figure('Name', 'Predicted vs Actual Severity', 'Position', [100 100 650 500]);
            scatter(R.yTest, R.yPred, 45, 'filled', 'MarkerFaceColor', [0.12 0.47 0.71], ...
                'MarkerEdgeColor', 'k', 'MarkerFaceAlpha', 0.7);
            hold on; lo = min([R.yTest; R.yPred]); hi = max([R.yTest; R.yPred]);
            plot([lo hi], [lo hi], 'r--', 'LineWidth', 2); grid on; axis square;
            ResultPlotter.setLabels('True Fault Severity (%)', 'Predicted Fault Severity (%)', ...
                'SVR Severity Estimation: Predicted vs. Actual');
            annotation('textbox', [0.16 0.72 0.28 0.15], 'String', {sprintf('R^2 = %.4f', R.r2), ...
                sprintf('RMSE = %.4f %%', R.rmse), sprintf('MAE = %.4f %%', R.mae)}, ...
                'FitBoxToText', 'on', 'BackgroundColor', [0.95 0.95 0.95], ...
                'EdgeColor', [0.7 0.7 0.7], 'FontSize', 10);
        end

        function residuals(R)
            figure('Name', 'Residual Analysis', 'Position', [150 150 1100 450]);
            subplot(1, 2, 1);
            scatter(R.yPred, R.err, 40, 'filled', 'MarkerFaceColor', [0.85 0.33 0.10], ...
                'MarkerEdgeColor', 'k', 'MarkerFaceAlpha', 0.6);
            hold on; yline(0, 'r--', 'LineWidth', 1.8); grid on;
            ResultPlotter.setLabels('Predicted Severity (%)', 'Residual (True - Predicted)', 'Residuals vs. Predicted Severity');
            subplot(1, 2, 2);
            histogram(R.err, 20, 'Normalization', 'pdf', 'FaceColor', [0.47 0.67 0.19], ...
                'EdgeColor', 'w', 'FaceAlpha', 0.75);
            hold on; g = linspace(min(R.err), max(R.err), 200);
            plot(g, normpdf(g, mean(R.err), std(R.err)), 'k-', 'LineWidth', 2); grid on;
            ResultPlotter.setLabels('Residual Error', 'Probability Density', 'Residual Distribution');
            legend({'Residuals', 'Gaussian Fit'}, 'Location', 'northeast');
        end

        function orderedTracking(R)
            figure('Name', 'Ordered Sample Tracking', 'Position', [200 200 1000 450]);
            [ys, si] = sort(R.yTest); yp = R.yPred(si); n = numel(ys);
            plot(1:n, ys, 'b-o', 'LineWidth', 1.5, 'MarkerSize', 4, 'MarkerFaceColor', 'b'); hold on;
            plot(1:n, yp, 'r--s', 'LineWidth', 1.5, 'MarkerSize', 4, 'MarkerFaceColor', 'r'); grid on;
            ResultPlotter.setLabels('Test Sample Index (Sorted by True Severity)', 'Fault Severity (%)', ...
                'True vs. Predicted Severity Across Test Samples');
            legend({'True Severity', 'SVR Prediction'}, 'Location', 'northwest');
        end

        function cvFolds(R)
            figure('Name', 'CV Fold Performance', 'Position', [250 250 750 450]);
            [fr, fm] = OutOfFoldPredictor.foldErrors(R.XtrS, R.yTrain, R.parts, R.params, R.weighter);
            b = bar([fr', fm'], 'grouped'); b(1).FaceColor = [0.2 0.4 0.8]; b(2).FaceColor = [0.8 0.4 0.2]; hold on;
            yline(mean(fr), 'b--', sprintf('Mean RMSE: %.3f', mean(fr)), 'LineWidth', 1.5);
            yline(mean(fm), 'r--', sprintf('Mean MAE: %.3f', mean(fm)), 'LineWidth', 1.5); grid on;
            ResultPlotter.setLabels('CV Fold Iteration Index', 'Error (%)', ...
                sprintf('Cross-Validation Performance Across %d Folds (%d Repeats)', ...
                R.cfg.FINAL_SEARCH_FOLDS, R.cfg.FINAL_SEARCH_REPEATS));
            legend({'RMSE', 'MAE'}, 'Location', 'northeast');
        end

        function generalizationGap(R)
            figure('Name', 'Overfitting Gap', 'Position', [300 300 650 450]);
            b = bar([R.trainRmse, R.rmse; R.trainR2, R.r2], 'grouped');
            b(1).FaceColor = [0.2 0.6 0.4]; b(2).FaceColor = [0.7 0.2 0.2];
            set(gca, 'XTickLabel', {'RMSE (%)', 'R^2 Score'}); grid on;
            ylabel('Metric Value', 'FontSize', 11, 'FontWeight', 'bold');
            title('Training vs. Test Performance Comparison', 'FontSize', 12);
            legend({'In-Sample Train', 'Out-of-Sample Test'}, 'Location', 'northeast');
            text(1, max([R.trainRmse, R.rmse]) * 0.8, sprintf('Gap: %.4f', R.rmse - R.trainRmse), ...
                'HorizontalAlignment', 'center', 'FontWeight', 'bold', 'FontSize', 10);
            text(2, min([R.trainR2, R.r2]) * 0.8, sprintf('Gap: %.4f', R.trainR2 - R.r2), ...
                'HorizontalAlignment', 'center', 'FontWeight', 'bold', 'FontSize', 10);
        end

        function permutationImportance(R)
            figure('Name', 'Feature Importance', 'Position', [350 350 850 550]);
            nf = size(R.XteS, 2); imp = zeros(1, nf);
            for f = 1:nf
                Xp = R.XteS; pr = zeros(1, 5);
                for rep = 1:5
                    Xp(:, f) = Xp(randperm(size(Xp, 1)), f);
                    pr(rep) = sqrt(mean((R.yTest - SvrTrainer.estimate(R.model, Xp, true)).^2));
                end
                imp(f) = mean(pr) - R.rmse;
            end
            [si, ii] = sort(imp, 'descend'); k = min(15, nf);
            barh(k:-1:1, si(1:k), 'FaceColor', [0.3 0.5 0.7]);
            yticks(1:k); yticklabels(R.names(ii(k:-1:1)));
            xlabel('Increase in Test RMSE upon Permutation', 'FontSize', 11, 'FontWeight', 'bold');
            title(sprintf('Top %d Permutation Feature Importances', k), 'FontSize', 12); grid on;
        end

        function searchHistory(R)
            figure('Name', 'Hyperparameter Search History', 'Position', [400 400 750 450]);
            rm = [R.trialLog.rmse]; feas = [R.trialLog.sv_ratio] <= R.cfg.TARGET_SV_RATIO;
            tgt = R.cfg.TARGET_SV_RATIO * 100;
            plot(find(feas), rm(feas), 'go', 'MarkerSize', 7, 'MarkerFaceColor', 'g', ...
                'DisplayName', sprintf('Feasible (SV \\le %.0f%%)', tgt)); hold on;
            plot(find(~feas), rm(~feas), 'rx', 'MarkerSize', 7, 'LineWidth', 1.5, ...
                'DisplayName', sprintf('Infeasible (SV > %.0f%%)', tgt)); grid on;
            ResultPlotter.setLabels('Search Evaluation Index', 'Validation RMSE (%)', 'SVR Hyperparameter Search History');
            legend('Location', 'northeast');
        end

        function multiclassDiagnostics(R)
            figure('Name', 'Multiclass Diagnostics', 'Position', [450 450 1100 500]);
            subplot(1, 2, 1);
            trc = histcounts(discretize(R.yTrain, R.edges), 1:(R.K+1));
            tec = histcounts(discretize(R.yTest,  R.edges), 1:(R.K+1));
            b = bar([trc; tec]', 'grouped'); b(1).FaceColor = [0.2 0.3 0.6]; b(2).FaceColor = [0.8 0.5 0.2];
            set(gca, 'XTickLabel', R.classNames); xtickangle(25);
            ylabel('Sample Count', 'FontSize', 11, 'FontWeight', 'bold');
            title('Dataset Class Distribution Across Splits', 'FontSize', 12);
            legend({'Train Split', 'Test Split'}, 'Location', 'northeast'); grid on;
            subplot(1, 2, 2);
            confusionchart(R.Cf, R.classNames, 'Title', sprintf('%d-Class Severity Confusion Matrix', R.K), ...
                'ColumnSummary', 'column-normalized', 'RowSummary', 'row-normalized');
        end

        function perClassMetrics(R)
            figure('Name', 'Per-Class Metrics', 'Position', [500 500 850 450]);
            mf = R.mf; b = bar([mf.prec(:), mf.rec(:), mf.spec(:), mf.f1(:)] * 100, 'grouped');
            cols = [0.12 0.47 0.71; 1 0.5 0.05; 0.17 0.63 0.17; 0.84 0.15 0.16];
            for i = 1:4, b(i).FaceColor = cols(i, :); end
            set(gca, 'XTickLabel', R.classNames); xtickangle(20);
            ylabel('Metric Value (%)', 'FontSize', 11, 'FontWeight', 'bold');
            title('Per-Class Diagnostic Performance Metrics (Test Set)', 'FontSize', 12);
            legend({'Precision', 'Recall/Sensitivity', 'Specificity', 'F1-Score'}, 'Location', 'southeast');
            grid on; ylim([0 105]);
        end

        function rocCurves(R)
            figure('Name', 'Multiclass ROC Curves', 'Position', [550 550 700 500]);
            cols = lines(R.K); hold on;
            for c = 1:R.K
                bt = (R.trueClass == c);
                sp = -abs(R.yPred - (R.edges(c) + R.edges(c+1)) / 2);
                sp(isinf(sp)) = R.yPred(isinf(sp));
                [xr, yr, ~, auc] = perfcurve(bt, sp, true);
                plot(xr, yr, 'LineWidth', 2, 'Color', cols(c, :), ...
                    'DisplayName', sprintf('%s (AUC = %.3f)', R.classNames{c}, auc));
            end
            plot([0 1], [0 1], 'k--', 'LineWidth', 1.5, 'DisplayName', 'Chance (AUC = 0.500)');
            grid on; axis square;
            ResultPlotter.setLabels('False Positive Rate (1 - Specificity)', 'True Positive Rate (Sensitivity)', ...
                'One-vs-Rest ROC Curves for Severity Classes');
            legend('Location', 'southeast');
            fprintf('All requested diagnostic and evaluation figures generated successfully.\n');
        end
    end
end
