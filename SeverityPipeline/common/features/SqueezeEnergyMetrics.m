classdef SqueezeEnergyMetrics
    % Cyclic (windowed) reconstructed energy, wavelet entropy and sub-band energy ratio
    methods (Static)
        function [s, e] = windowRanges(t, windowSec)
            dt = t(2) - t(1); N = length(t);
            spw = round(windowSec / dt);
            c = (1:max(1, ceil(N / spw)))';
            s = (c - 1) * spw + 1;
            e = min(c * spw, N);
        end
        function E = cyclicEnergy(instE, t, windowSec)
            [s, e] = SqueezeEnergyMetrics.windowRanges(t, windowSec);
            E = zeros(numel(s), 1);
            for c = 1:numel(s), E(c) = sum(instE(s(c):e(c))); end
        end
        function H = waveletEntropy(Tx, t, windowSec)
            [s, e] = SqueezeEnergyMetrics.windowRanges(t, windowSec);
            H = zeros(numel(s), 1);
            for w = 1:numel(s)
                H(w) = SqueezeEnergyMetrics.shannon(sum(Tx(:, s(w):e(w)).^2, 2));
            end
        end
        function H = shannon(E)
            E = E(:); total = sum(E);
            if total <= 0, H = 0; return; end
            p = E / total; p(p == 0) = [];
            H = -sum(p .* log2(p));
        end
        function r = subbandEnergyRatio(Tx, bandIdx, lowerHist, t, windowSec)
            NumTime = size(Tx, 2);
            bandE = zeros(NumTime, 1); totE = zeros(NumTime, 1);
            for idx = 1:NumTime
                coefs = Tx(bandIdx{idx}, idx);
                coefs(coefs < lowerHist(idx)) = 0;
                bandE(idx) = sum(coefs.^2);
                totE(idx)  = sum(Tx(:, idx).^2);
            end
            [s, e] = SqueezeEnergyMetrics.windowRanges(t, windowSec);
            r = zeros(numel(s), 1);
            for w = 1:numel(s)
                b = sum(bandE(s(w):e(w))); tt = sum(totE(s(w):e(w)));
                if tt > 0, r(w) = b / tt; end
            end
        end
    end
end
