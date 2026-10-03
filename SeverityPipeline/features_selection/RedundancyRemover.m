classdef RedundancyRemover
    % Drops a group when its window-by-window |correlation| with a more relevant group is high
    methods (Static)
        function final = remove(X, groups, retained, rel, cfg, verbose)
            C = abs(corr(X));                         % scaling does not change correlations
            [~, loc] = ismember(retained, rel.group);
            [~, o] = sort(rel.relevance_score(loc), 'descend');
            ordered = retained(o);
            removed = {}; map = {};
            for i = 1:numel(ordered)
                gi = ordered{i};
                if any(strcmp(removed, gi)), continue; end
                for j = i+1:numel(ordered)
                    gj = ordered{j};
                    if any(strcmp(removed, gj)) || any(strcmp(cfg.ProtectedGroups, gj)), continue; end
                    if median(RedundancyRemover.pairCorr(C, groups, gi, gj)) >= cfg.RedundancyThreshold
                        removed{end+1} = gj; map(end+1, :) = {gj, gi}; %#ok<AGROW>
                    end
                end
            end
            final = retained(~ismember(retained, removed));
            if verbose
                SelectionRules.header('STEP 4 - REDUNDANCY REMOVAL');
                fprintf('Redundancy threshold: %g\n\nRemoved group -> retained representative\n%s\n', cfg.RedundancyThreshold, repmat('-', 1, 42));
                for k = 1:size(map, 1), fprintf('%s -> %s\n', map{k, 1}, map{k, 2}); end
                fprintf('\nBefore redundancy removal: %d groups\nAfter redundancy removal: %d groups\nRemoved: %d groups\n', ...
                    numel(retained), numel(final), numel(removed));
                SelectionRules.list('RETAINED NON-REDUNDANT GROUPS', final);
            end
        end
        function c = pairCorr(C, groups, a, b)
            ia = groups.Index(a); ib = groups.Index(b);
            [w, ka, kb] = intersect(groups.Wins{ia}, groups.Wins{ib});
            if ~isempty(w), c = C(sub2ind(size(C), groups.Cols{ia}(ka), groups.Cols{ib}(kb)));
            else,           c = reshape(C(groups.Cols{ia}, groups.Cols{ib}), [], 1);
            end
        end
    end
end