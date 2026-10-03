classdef NumpyStats
    % numpy-style linear-interpolation percentile (MATLAB's prctile uses another convention)
    methods (Static)
        function v = percentile(x, q)             % q in [0,100], scalar or vector
            xs = sort(x(:)); n = numel(xs);
            if n == 1, v = repmat(xs, size(q)); return; end
            v = interp1(1:n, xs, q/100*(n-1) + 1);
        end
    end
end