classdef StandardScaler
    % z-score scaling (fitted on train only)
    properties (SetAccess = private)
        Mu
        Sigma
    end
    methods
        function obj = fit(obj, X)
            obj.Mu = mean(X, 1);
            obj.Sigma = std(X, 0, 1);
            obj.Sigma(obj.Sigma == 0) = 1;
        end
        function Xo = transform(obj, X), Xo = (X - obj.Mu) ./ obj.Sigma; end
    end
end
