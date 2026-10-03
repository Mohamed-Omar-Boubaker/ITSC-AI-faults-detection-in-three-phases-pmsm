classdef SeverityExperiment
    % Orchestrator of the whole training workflow (same steps/order as the original script)
    properties (SetAccess = immutable)
        Cfg
    end
    methods
        function obj = SeverityExperiment(cfg), obj.Cfg = cfg; end

        function run(obj)
            cfg = obj.Cfg; hb = cfg.HEALTHY_THRESHOLD;
            rng(42, 'twister'); maxNumCompThreads(1);
            weighter = SampleWeighter(cfg);

            %% 1-2. Load + stratified 70/30 split
            [X, y, rawNames] = SeverityDataset.fromMat(cfg.DATA_PATH);
            fprintf('Loaded %d samples, %d raw features, %d unique severity levels\n', numel(y), size(X, 2), numel(unique(y)));
            cv = StratifiedSplitter.holdout(y, 0.3);
            yTr = y(cv.training); yTe = y(cv.test);
            fprintf('\nTrain samples: %d, Test samples: %d\n', numel(yTr), numel(yTe));

            %% 3 + 5. Preprocessing fitted on TRAIN only (family collapse, zero-variance filter)
            collapser = FeatureFamilyCollapser().fit(X(cv.training, :), rawNames);
            XtrC = collapser.transform(X(cv.training, :), rawNames);
            XteC = collapser.transform(X(cv.test, :), rawNames);
            fprintf('Collapsed to %d features from %d families\n', size(XtrC, 2), numel(collapser.FamilyRoots));
            vf = VarianceFilter().fit(XtrC);
            Xtr = vf.transform(XtrC); Xte = vf.transform(XteC);
            names = collapser.OutputNames(vf.Valid);
            fprintf('Healthy samples in train: %d/%d\n', sum(yTr == 0), numel(yTr));

            %% 6. Optional nested CV
            if cfg.RUN_NESTED_CV, NestedCrossValidator(cfg).run(Xtr, yTr); end

            %% 7. Final hyper-parameter search (SV-ratio constrained)
            fprintf('\n================ FINAL HYPERPARAMETER SEARCH ================\n');
            wTr = weighter.compute(yTr);
            scaler = StandardScaler().fit(Xtr);
            XtrS = scaler.transform(Xtr); XteS = scaler.transform(Xte);
            parts = CrossValidation.repeatedKFold(numel(yTr), cfg.FINAL_SEARCH_FOLDS, cfg.FINAL_SEARCH_REPEATS);
            [best, trialLog] = ConstrainedSvrSearch(cfg).run(XtrS, yTr, wTr, parts, cfg.FINAL_SEARCH_MAX_EVALS);
            ConsoleReport.hyperparameters(best, trialLog, cfg);

            %% 8. Out-of-fold predictions
            [oof, cnt] = OutOfFoldPredictor.compute(XtrS, yTr, parts, best, weighter);
            fprintf('\nOOF predictions generated for %d/%d training samples (avg %.1f repeats/sample)\n', sum(cnt > 0), numel(yTr), mean(cnt));
            sens = sum((oof > hb) & (yTr > 0)) / max(sum(yTr > 0), 1);
            spec = sum(~(oof > hb) & (yTr == 0)) / max(sum(yTr == 0), 1);
            fprintf('Healthy Cutoff Boundary: Strictly Fixed at %.2f%%\n', hb);
            fprintf('OOF Health Detection Performance: Sensitivity = %.2f%%, Specificity = %.2f%%\n', sens*100, spec*100);

            %% 9. Final model
            model   = SvrTrainer.train(XtrS, yTr, best, wTr);
            yPredTr = SvrTrainer.estimate(model, XtrS, true);
            yPred   = SvrTrainer.estimate(model, XteS, true);
            nSv = size(model.SupportVectors, 1); svRatio = nSv / numel(yTr);

            %% 10-12. Metrics
            m = RegressionMetrics.compute(yTe, yPred);
            ConsoleReport.regression(m, numel(yTe));
            ConsoleReport.severity(RegressionMetrics.severityBreakdown(yTe, yPred, cfg.SEVERITY_THRESHOLD));
            bin = BinaryHealthMetrics.compute(yTe, yPred, hb);
            ConsoleReport.binary(bin, hb);
            mt = RegressionMetrics.compute(yTr, yPredTr);
            ConsoleReport.diagnostics(mt.rmse, m.rmse, mt.r2, m.r2, nSv, numel(yTr));

            %% 13. Sanity checks, SV analysis, noise robustness, learning curve, calibration
            ModelDiagnostics(cfg, model, best, weighter, XtrS, yTr, XteS, yTe, y).runAll(m, yPred);

            %% 14. Multiclass boundaries (DP) -> frozen -> evaluated once on the test set
            fprintf('\n================ MULTICLASS BOUNDARY OPTIMIZATION ================\n');
            ev  = BoundaryEvaluator(hb, cfg.SCORE_WEIGHTS, cfg.USE_BALANCED_ACC);
            opt = BoundaryOptimizerDP(ev, cfg.TARGET_CLASS_WIDTH, cfg.WIDTH_BONUS_WEIGHT, cfg.N_CANDIDATE_GRID);
            sel = BoundarySelector(cfg, opt, BoundaryCrossValidator(cfg, opt, ev)).select(yTr, oof, bin);
            K = sel.NClasses;
            testRes = ev.evaluate(yTe, yPred, sel.Boundaries, K, 0, Inf);
            ConsoleReport.multiclassTest(testRes, sel);

            %% 14c. Confusion matrices (raw levels, binary, chosen classes)
            fprintf('\n================ CONFUSION MATRICES ================\n');
            levels = sort(unique(y));
            lvlLabels = arrayfun(@(v) sprintf('%.2f', v), levels, 'UniformOutput', false);
            [~, predLvl] = min(abs(yPred - levels'), [], 2);
            [~, trueLvl] = min(abs(yTe   - levels'), [], 2);
            [C22, m22] = ConfusionPlotter.suite(trueLvl, predLvl, lvlLabels, 'Confusion Matrix - 22 Severity Levels');
            ConsoleReport.perClass('Per-Level Test Metrics (22 severity levels)', lvlLabels, m22, C22, 10);
            binLabels = {sprintf('Healthy (0-%.2f%%)', hb), sprintf('Faulty (>%.2f%%)', hb)};
            ConfusionPlotter.suite(double(bin.actualFaulty) + 1, double(bin.predictedFaulty) + 1, ...
                binLabels, 'Confusion Matrix - Binary (Healthy vs Faulty)');
            trueClass = discretize(yTe, sel.Edges); predClass = discretize(yPred, sel.Edges);
            ConfusionPlotter.suite(trueClass, predClass, sel.ClassNames, sprintf('Confusion Matrix - %d Chosen Classes', K));

            %% 15-16. Save model bundle + target summary
            ModelBundleWriter.write(cfg.SAVE_PATH, struct('model', model, 'collapser', collapser, ...
                'varianceFilter', vf, 'scaler', scaler, 'rawFeatureNames', {rawNames}, 'featureNames', {names}, ...
                'bestParams', best, 'cfg', cfg.toStruct(), 'edges', sel.Edges, 'classNames', {sel.ClassNames}));
            ConsoleReport.targets(cfg, svRatio, bin, testRes.accuracy, K);

            %% 17. Figures
            fprintf('\n==================================================\n GENERATING PERFORMANCE & DIAGNOSTIC PLOTS\n==================================================\n');
            R = struct('yTest', yTe, 'yPred', yPred, 'yTrain', yTr, 'err', m.err, 'rmse', m.rmse, 'mae', m.mae, 'r2', m.r2, ...
                'trainRmse', mt.rmse, 'trainR2', mt.r2, 'cfg', cfg, 'parts', {parts}, 'XtrS', XtrS, 'XteS', XteS, ...
                'model', model, 'params', best, 'weighter', weighter, 'trialLog', trialLog, 'edges', sel.Edges, ...
                'classNames', {sel.ClassNames}, 'K', K, 'Cf', testRes.C, 'mf', testRes.metrics, ...
                'trueClass', trueClass, 'names', {names});
            ResultPlotter.all(R);
        end
    end
end
