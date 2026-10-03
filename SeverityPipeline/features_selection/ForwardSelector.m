classdef ForwardSelector
    % Sequential forward selection of whole groups (protected groups seeded first)
    properties (SetAccess = immutable)
        Groups
        Evaluator
        Cfg
    end
    methods
        function obj = ForwardSelector(groups, evaluator, cfg)
            obj.Groups = groups; obj.Evaluator = evaluator; obj.Cfg = cfg;
        end
        function hist = select(obj, X, y, candidates, rel, maxGroups, verbose)
            cfg = obj.Cfg; G = obj.Groups;
            if verbose, SelectionRules.header('STEP 5 - SEQUENTIAL FORWARD GROUP SELECTION'); end
            prot = cfg.ProtectedGroups(ismember(cfg.ProtectedGroups, candidates));
            remaining = candidates(~ismember(candidates, prot));
            selected = {}; rows = {}; it = 1; bestRmse = Inf;
            if ~isempty(prot)
                selected = prot;
                r0 = obj.Evaluator.cvRmse(X, y, G.cols(selected));
                for k = 1:numel(prot), rows(end+1, :) = {0, prot{k}, numel(selected), r0}; end %#ok<AGROW>
                if verbose, fprintf('\nSeeded protected groups: %s\nCV RMSE with protected groups only: %.4f\n', strjoin(prot, ', '), r0); end
                bestRmse = r0;
            end
            if isempty(maxGroups), lim = numel(candidates); else, lim = maxGroups; end
            lim = max(lim, numel(prot)); noImp = 0;
            while ~isempty(remaining) && numel(selected) < lim
                pool = remaining;
                if isfinite(cfg.CandidatesPerIter)             % test only the top-K by relevance
                    [~, loc] = ismember(pool, rel.group);
                    [~, o] = sort(rel.relevance_score(loc), 'descend');
                    pool = pool(o(1:min(cfg.CandidatesPerIter, numel(o))));
                end
                bestG = ''; bestR = Inf;
                for k = 1:numel(pool)
                    r = obj.Evaluator.cvRmse(X, y, G.cols([selected, pool(k)]));
                    if r < bestR, bestR = r; bestG = pool{k}; end
                end
                if isempty(bestG), bestG = pool{1}; end
                selected{end+1} = bestG; remaining(strcmp(remaining, bestG)) = []; %#ok<AGROW>
                rows(end+1, :) = {it, bestG, numel(selected), bestR}; %#ok<AGROW>
                if verbose
                    fprintf('\nIteration %d (tested %d candidates)\nBest group added: %s\nCV RMSE: %.4f\nSelected groups: %d | columns: %d\n%s\n', ...
                        it, numel(pool), bestG, bestR, numel(selected), numel(G.cols(selected)), repmat('-', 1, 62));
                end
                if bestR < bestRmse - 1e-4, bestRmse = bestR; noImp = 0; else, noImp = noImp + 1; end
                if noImp >= cfg.Patience
                    if verbose, fprintf('No improvement for %d iterations, stopping early\n', cfg.Patience); end
                    break;
                end
                it = it + 1;
            end
            hist = cell2table(rows, 'VariableNames', {'iteration', 'group_added', 'n_groups', 'cv_rmse'});
        end
    end
end