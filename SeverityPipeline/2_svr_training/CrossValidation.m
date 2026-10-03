classdef CrossValidation
    methods (Static)
        function parts = repeatedKFold(n, K, R)
            parts = cell(R, 1);
            for r = 1:R, parts{r} = cvpartition(n, 'KFold', K); end
        end
    end
end
