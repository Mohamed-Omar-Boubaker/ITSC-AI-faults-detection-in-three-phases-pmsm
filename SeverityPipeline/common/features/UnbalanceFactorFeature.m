classdef UnbalanceFactorFeature < WindowFeature
    % Current unbalance factor: negative / positive sequence RMS ratio (%)
    methods
        function fs = extract(~, w)
            v = (sqrt(mean(abs(w.INeg(:)).^2)) / (sqrt(mean(abs(w.IPos(:)).^2)) + eps)) * 100;
            fs = FeatureSet({'Cross_CUF'}, v);
        end
    end
end
