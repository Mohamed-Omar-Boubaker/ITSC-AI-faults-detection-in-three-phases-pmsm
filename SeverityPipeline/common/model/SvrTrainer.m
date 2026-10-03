classdef SvrTrainer
    % Gaussian-kernel epsilon-SVR wrapper (shared by training and inference)
    methods (Static)
        function mdl = train(X, y, p, w)
            mdl = fitrsvm(X, y, 'KernelFunction', 'gaussian', 'BoxConstraint', p.BoxConstraint, ...
                'KernelScale', p.KernelScale, 'Epsilon', p.Epsilon, 'Weights', w, 'Standardize', false);
        end
        function yhat = estimate(mdl, X, clipAtZero)
            yhat = predict(mdl, X);
            if clipAtZero, yhat = max(0, yhat); end
        end
    end
end
