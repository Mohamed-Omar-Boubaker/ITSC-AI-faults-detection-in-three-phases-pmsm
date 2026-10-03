classdef VarianceFilter
    % Drops zero-variance columns (fitted on train only)
    properties (SetAccess = private)
        Valid
    end
    methods
        function obj = fit(obj, X), obj.Valid = std(X, 0, 1) > 0; end
        function Xo = transform(obj, X), Xo = X(:, obj.Valid); end
    end
end
