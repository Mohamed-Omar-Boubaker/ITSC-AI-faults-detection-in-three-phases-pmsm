classdef SelectedModel
    % Fitted model + the columns (and optional scaler) it needs
    properties (SetAccess = immutable)
        Model
        Cols
        Mu
        Sd
    end
    methods
        function obj = SelectedModel(model, cols, mu, sd)
            obj.Model = model; obj.Cols = cols; obj.Mu = mu; obj.Sd = sd;
        end
        function yhat = estimate(obj, X)
            Z = X(:, obj.Cols);
            if ~isempty(obj.Mu), Z = (Z - obj.Mu) ./ obj.Sd; end
            yhat = predict(obj.Model, Z);
        end
    end
end