classdef SelectionRules
    % Threshold rules: relevance filtering and final subset size, plus console helpers
    methods (Static)
        function retained = relevance(rel, cfg, scorer, verbose)
            thr = NumpyStats.percentile(rel.relevance_score, 100 - cfg.RelevancePercentile);
            retained = rel.group(rel.relevance_score >= thr)';
            if numel(retained) < cfg.MinRetainedGroups
                retained = rel.group(1:min(cfg.MinRetainedGroups, height(rel)))';
            end
            forced = {};
            for k = 1:numel(cfg.ProtectedGroups)
                g = cfg.ProtectedGroups{k};
                if any(strcmp(rel.group, g)) && ~any(strcmp(retained, g))
                    retained{end+1} = g; forced{end+1} = g; %#ok<AGROW>
                end
            end
            removed = setdiff(rel.group', retained, 'stable');
            if verbose
                SelectionRules.header('STEP 3 - GROUP RELEVANCE FILTERING');
                show = rel(:, scorer.displayVars()); show.Rank = (1:height(show))';
                disp(movevars(show, 'Rank', 'Before', 1));
                fprintf('\nRELEVANCE FILTER RESULT\nOriginal groups: %d\nRetained groups: %d\nRemoved groups: %d\n', ...
                    height(rel), numel(retained), numel(removed));
                if ~isempty(forced), fprintf('Force-retained protected groups (below threshold): %s\n', strjoin(forced, ', ')); end
                SelectionRules.list('RETAINED BASE FEATURES', retained);
                SelectionRules.list('REMOVED BASE FEATURES', removed);
            end
        end

        function n = subsetSize(hist, protectedCount, cfg)
            SelectionRules.header('FINAL SUBSET SIZE SELECTION');
            minR = min(hist.cv_rmse); thr = minR * (1 + cfg.RmseTolerance);
            elig = sortrows(hist(hist.cv_rmse <= thr, :), 'n_groups');
            chosen = elig.n_groups(1);
            n = min(max([chosen, protectedCount, cfg.MinFinalGroups]), max(hist.n_groups));
            fprintf('Best CV RMSE: %.4f\n%.0f%% tolerance threshold: %.4f\n', minR, cfg.RmseTolerance*100, thr);
            fprintf('Tolerance-based size: %d | Protected floor: %d | Min floor: %d\n', chosen, protectedCount, cfg.MinFinalGroups);
            fprintf('\nSelected subset:\n%d base feature groups\n', n);
            fprintf('\nReason:\nmax(smallest size within RMSE tolerance, protected group count, minimum floor)\n');
        end

        function header(t), fprintf('%s\n%s\n%s\n', repmat('=', 1, 66), t, repmat('=', 1, 66)); end
        function list(title, items)
            fprintf('\n%s:\n', title);
            for i = 1:numel(items), fprintf('%d. %s\n', i, items{i}); end
        end
    end
end