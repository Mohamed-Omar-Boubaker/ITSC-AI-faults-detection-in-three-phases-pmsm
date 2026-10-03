classdef SvrFactory
    % sklearn-equivalent SVR pieces (rbf kernel, population-std scaler, gamma='scale')
    methods (Static)
        function [mu, sd] = scalerFit(X), mu = mean(X, 1); sd = std(X, 1, 1); sd(sd == 0) = 1; end
        function g = gammaScale(Xs), g = 1 / (size(Xs, 2) * var(Xs(:), 1)); end
        function m = fit(Xs, y, C, epsilon, gamma)
            m = fitrsvm(Xs, y, 'KernelFunction', 'gaussian', 'KernelScale', 1/sqrt(gamma), ...
                'BoxConstraint', C, 'Epsilon', epsilon, 'Standardize', false);
        end
    end
end