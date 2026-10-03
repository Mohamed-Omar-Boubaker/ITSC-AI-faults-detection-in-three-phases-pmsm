classdef ForestRelevance < RelevanceScorer
    % relevance = normalised median random-forest importance of the group's windows
    properties (SetAccess = immutable)
        Trees
        Depth
    end
    methods
        function obj = ForestRelevance(trees, depth), obj.Trees = trees; obj.Depth = depth; end
        function v = displayVars(~), v = {'group', 'imp_mean', 'imp_median', 'imp_max', 'relevance_score'}; end
        function rel = score(obj, X, y, groups)
            m = ForestFactory.fit(X, y, obj.Trees, obj.Depth, 2, 1, 'all');
            imp = predictorImportance(m);
            G = groups.numGroups(); a = zeros(G, 1); b = a; c = a;
            for g = 1:G
                v = imp(groups.Cols{g}); a(g) = mean(v); b(g) = median(v); c(g) = max(v);
            end
            r = (b - min(b)) / (max(b) - min(b) + 1e-12);
            rel = table(groups.Bases(:), a, b, c, r, ...
                'VariableNames', {'group', 'imp_mean', 'imp_median', 'imp_max', 'relevance_score'});
            rel = sortrows(rel, 'relevance_score', 'descend');
        end
    end
end