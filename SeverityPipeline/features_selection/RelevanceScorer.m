classdef (Abstract) RelevanceScorer
    % Scores every feature group; returns a table sorted by 'relevance_score' (descending)
    methods (Abstract)
        rel = score(obj, X, y, groups)
        v = displayVars(obj)
    end
end