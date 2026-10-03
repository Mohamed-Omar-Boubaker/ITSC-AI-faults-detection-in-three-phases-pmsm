classdef FeatureSelectionPipeline
    % Orchestrator: relevance -> redundancy -> forward selection -> size -> stability -> final model -> save
    properties (SetAccess = immutable)
        Cfg
        Scorer       % RelevanceScorer
        Evaluator    % SubsetEvaluator
        Tuner        % object with fit(X, y, cols) -> [SelectedModel, params]
    end
    methods
        function obj = FeatureSelectionPipeline(cfg, scorer, evaluator, tuner)
            obj.Cfg = cfg; obj.Scorer = scorer; obj.Evaluator = evaluator; obj.Tuner = tuner;
        end

        function run(obj)
            cfg = obj.Cfg; rng(cfg.Seed, 'twister');
            SelectionRules.header('STEP 1 - DATASET AND TRAIN/TEST SPLIT');
            S = load(cfg.InFile);
            featRaw = double(S.feat_table); names = cellstr(S.names(:))';

            SelectionRules.header('STEP 2 - FEATURE GROUP VALIDATION');
            groups = FeatureGroups(names, size(featRaw, 2), cfg.RequiredWindows);
            fprintf('\nDetected base feature groups: %d\nTotal columns: %d\n\nGroup validation: PASSED\n\nExample groups and window counts:\n', ...
                groups.numGroups(), size(featRaw, 2));
            for i = 1:min(5, groups.numGroups())
                fprintf('%d. %s -> %d windows\n', i, groups.Bases{i}, numel(groups.Cols{i}));
            end
            missing = cfg.ProtectedGroups(~cellfun(@(g) groups.has(g), cfg.ProtectedGroups));
            if ~isempty(missing), fprintf('\nWARNING: protected groups not found in dataset: %s\n', strjoin(missing, ', ')); end

            [X, y] = DatasetPreparer.clean(featRaw, S.percent);
            fprintf('\nOriginal samples: %d\nTotal columns: %d\nBase feature groups: %d\n', size(X, 1), size(X, 2), groups.numGroups());
            [Xtr, Xte, ytr, yte] = DatasetPreparer.splitAndImpute(X, y, cfg.Seed);
            fprintf('\nTraining samples: %d\nHeld-out test samples: %d\n', size(Xtr, 1), size(Xte, 1));
            fprintf('\nNumber of severity levels in training: %d\nNumber of severity levels in test: %d\n', numel(unique(ytr)), numel(unique(yte)));

            rel      = obj.Scorer.score(Xtr, ytr, groups);
            retained = SelectionRules.relevance(rel, cfg, obj.Scorer, true);
            nonRed   = RedundancyRemover.remove(Xtr, groups, retained, rel, cfg, true);
            hist     = ForwardSelector(groups, obj.Evaluator, cfg).select(Xtr, ytr, nonRed, rel, [], true);

            nProt  = numel(intersect(cfg.ProtectedGroups, nonRed, 'stable'));
            nFinal = SelectionRules.subsetSize(hist, nProt, cfg);
            ordered = sortrows(hist, 'n_groups');                   % stable: protected first, then addition order
            finalGroups = ordered.group_added(1:nFinal)';
            StabilityAnalyzer.run(Xtr, ytr, groups, cfg, obj.Scorer, obj.Evaluator, nonRed, rel, nFinal);

            [~, o] = sort(cellfun(@(g) groups.firstCol(g), finalGroups));
            finalOrdered = finalGroups(o);                          % original column order
            finalCols = groups.cols(finalOrdered);

            SelectionRules.header('STEP 7 - FINAL MODEL');
            fprintf('Final base feature groups: %d\nFinal feature columns: %d\n', numel(finalOrdered), numel(finalCols));
            SelectionRules.list('Selected features', finalOrdered);
            [model, params] = obj.Tuner.fit(Xtr, ytr, finalCols);
            fprintf('\nBest parameters:\n'); disp(params);

            SelectionRules.header('STEP 8 - FINAL HELD-OUT TEST PERFORMANCE');
            pr = model.estimate(Xte);
            fprintf('Test RMSE: %.4f\nTest MAE: %.4f\nTest R2: %.4f\n', sqrt(mean((yte - pr).^2)), mean(abs(yte - pr)), ...
                1 - sum((yte - pr).^2) / sum((yte - mean(yte)).^2));

            SelectionRules.header('FINAL FEATURE SELECTION REPORT');
            nCols = @(g) numel(groups.cols(g));
            fprintf('\nOriginal:\n%d base features\n%d columns\n', groups.numGroups(), size(X, 2));
            fprintf('\nAfter relevance filtering:\n%d base features\n%d columns\n', numel(retained), nCols(retained));
            fprintf('\nAfter redundancy removal:\n%d base features\n%d columns\n', numel(nonRed), nCols(nonRed));
            fprintf('\nFinal selected pool:\n%d base features\n%d columns\n', numel(finalOrdered), numel(finalCols));
            fprintf('\nColumn Reduction:\n%.1f%%\n', 100 * (1 - numel(finalCols) / size(X, 2)));
            SelectionRules.header('ALL WINDOWS OF EACH SELECTED FEATURE');
            for k = 1:numel(finalOrdered)
                g = finalOrdered{k}; w = groups.Wins{groups.Index(g)};
                fprintf('\nFeature: %s (%d windows)\n', g, numel(w));
                for i = 1:numel(w), fprintf('  %s_w%d\n', g, w(i)); end
            end
            SelectionRules.list('FINAL REMOVED FEATURES', setdiff(groups.Bases, finalOrdered, 'stable'));

            if cfg.SaveCleaned, FilteredDatabaseWriter.write(cfg.OutFile, X, names, finalCols, y, []);
            else,               FilteredDatabaseWriter.write(cfg.OutFile, featRaw, names, finalCols, S.percent, S.degree);
            end
        end
    end
end