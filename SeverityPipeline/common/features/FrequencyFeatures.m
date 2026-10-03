classdef FrequencyFeatures < FeatureGroup
    % 7 spectral features on the two-sided complex FFT magnitude
    properties (SetAccess = immutable)
        Fs
        F0
        H
    end
    methods
        function obj = FrequencyFeatures(fs, f0, H)
            obj.Fs = fs; obj.F0 = f0; obj.H = H;
        end
        function r = registry(obj)
            fs = obj.Fs; f0 = obj.F0; H = obj.H;
            r = { ...
                'Spectral_Centroid',   @(x) FrequencyFeatures.centroid(x, fs); ...
                'Spectral_Spread',     @(x) FrequencyFeatures.spread(x, fs); ...
                'Spectral_Kurtosis',   @(x) FrequencyFeatures.spectralKurtosis(x, fs); ...
                'Target_Band_Energy',  @(x) FrequencyFeatures.bandEnergy(x, fs, 2*f0-10, 2*f0+10); ...
                'Harmonic_Distortion', @(x) FrequencyFeatures.harmonicDistortion(x, fs, f0, H); ...
                'Spectral_Flatness',   @(x) FrequencyFeatures.flatness(x, fs); ...
                'Kurtogram_Peak',      @(x) FrequencyFeatures.kurtogramPeak(x, fs)};
        end
    end
    methods (Static)
        function [mag, f] = spectrum(x, fs)
            x = x(:); N = length(x);
            mag = abs(fftshift(fft(x)) / N);
            f = (-floor(N/2) : ceil(N/2)-1)' * (fs / N);
        end
        function v = centroid(x, fs)
            [m, f] = FrequencyFeatures.spectrum(x, fs);
            v = sum(f .* m) / (sum(m) + eps);
        end
        function v = spread(x, fs)
            [m, f] = FrequencyFeatures.spectrum(x, fs);
            fsc = FrequencyFeatures.centroid(x, fs);
            v = sqrt(sum(((f - fsc).^2) .* m) / (sum(m) + eps));
        end
        function v = spectralKurtosis(x, fs)
            m = FrequencyFeatures.spectrum(x, fs);
            v = mean(((m - mean(m)) / (std(m, 1) + eps)).^4);
        end
        function v = bandEnergy(x, fs, fLo, fHi)
            [m, f] = FrequencyFeatures.spectrum(x, fs);
            idx = (abs(f) >= fLo) & (abs(f) <= fHi);
            v = sum(m(idx).^2) / (sum(m.^2) + eps);
        end
        function v = harmonicDistortion(x, fs, f0, H)
            [m, f] = FrequencyFeatures.spectrum(x, fs);
            [~, i0] = min(abs(f - f0)); mag0 = m(i0);
            hsq = 0;
            for h = 2:H
                [~, ip] = min(abs(f - h * f0));
                [~, in_] = min(abs(f + h * f0));
                hsq = hsq + m(ip)^2 + m(in_)^2;
            end
            v = sqrt(hsq) / (mag0 + eps);
        end
        function v = flatness(x, fs)
            m = FrequencyFeatures.spectrum(x, fs);
            v = exp(mean(log(m + 1e-12))) / (mean(m) + eps);
        end
        function v = kurtogramPeak(x, fs)
            K = kurtogram(abs(x(:)), fs); v = max(K(:));
        end
    end
end
