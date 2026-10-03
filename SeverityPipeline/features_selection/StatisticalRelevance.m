classdef StatisticalRelevance < RelevanceScorer
    % relevance = 0.5*norm(median |Spearman|) + 0.5*norm(median mutual information)
    properties (SetAccess = immutable)
        Seed
    end
    methods
        function obj = StatisticalRelevance(seed), obj.Seed = seed; end
        function v = displayVars(~), v = {'group', 'spearman_median', 'mi_median', 'relevance_score'}; end
        function rel = score(obj, X, y, groups)
            sp = abs(corr(tiedrank(X), tiedrank(y(:))))';     % Spearman = Pearson on ranks
            sp(isnan(sp)) = 0;
            mi = MutualInfo.regression(X, y, 3, obj.Seed);
            G = groups.numGroups(); spMed = zeros(G, 1); miMed = zeros(G, 1);
            for g = 1:G
                c = groups.Cols{g}; spMed(g) = median(sp(c)); miMed(g) = median(mi(c));
            end
            nrm = @(s) (s - min(s)) / (max(s) - min(s) + 1e-12);
            rel = table(groups.Bases(:), spMed, miMed, 0.5*nrm(spMed) + 0.5*nrm(miMed), ...
                'VariableNames', {'group', 'spearman_median', 'mi_median', 'relevance_score'});
            rel = sortrows(rel, 'relevance_score', 'descend');
        end
    end
end