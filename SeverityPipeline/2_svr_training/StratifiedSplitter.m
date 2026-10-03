classdef StratifiedSplitter
    % Stratified hold-out split
    methods (Static)
        function cv = holdout(y, frac)
            if numel(unique(y)) < 30
                cv = cvpartition(y, 'Holdout', frac);
            else
                cv = cvpartition(discretize(y, 10), 'Holdout', frac);
            end
        end
    end
end
