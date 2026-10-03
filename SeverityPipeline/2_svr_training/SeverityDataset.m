classdef SeverityDataset
    % Loads the (filtered) feature table + severity labels
    methods (Static)
        function [X, y, names] = fromMat(path)
            S = load(path);
            X = S.feat_table;
            y = S.percent(:);
            if isfield(S, 'names')
                names = S.names;
            else
                names = arrayfun(@(i) sprintf('Feature_%d', i), 1:size(X, 2), 'UniformOutput', false);
            end
        end
    end
end
