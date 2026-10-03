classdef MorseCwt
    % Bilateral (positive / negative frequency) Morse wavelet transform, FFT domain
    properties (SetAccess = immutable)
        Cfg
    end
    methods
        function obj = MorseCwt(cfg), obj.Cfg = cfg; end

        function [scales, f] = morseScales(obj, N)
            c = obj.Cfg;
            w0 = (c.BetaOpt / c.Gamma)^(1 / c.Gamma);
            fmax = c.Fs / 2; fmin = c.Fs / N;
            Nv = floor(log2(fmax / fmin) * c.VoicesPerOctave);
            f = fmax ./ 2.^((0:Nv) / c.VoicesPerOctave);
            f = f(f >= fmin).';
            scales = (w0 * c.Fs) ./ (2 * pi * f);
        end

        function [W, dW] = transform(obj, x, scales)
            c = obj.Cfg;
            x = x(:).'; N = length(x);
            X = fft(x);
            w = (2*pi/N) * (0:N-1);
            omega = (2*pi*c.Fs/N) * [0:floor((N-1)/2), -floor(N/2):-1];
            Ns = length(scales);
            W = zeros(Ns, N, 2); dW = zeros(Ns, N, 2);
            for k = 1:Ns
                a = scales(k);
                psiP = MorseCwt.waveletFFT(w, a, c.BetaOpt, c.Gamma, 'pos');
                W(k,:,1)  = ifft(X .* conj(psiP));
                dW(k,:,1) = ifft(X .* (1i * omega) .* conj(psiP));
                psiN = MorseCwt.waveletFFT(w, a, c.BetaOpt, c.Gamma, 'neg');
                W(k,:,2)  = ifft(X .* conj(psiN));
                dW(k,:,2) = ifft(X .* (1i * omega) .* conj(psiN));
            end
        end

        function C = admissibilityConstant(obj)
            b = obj.Cfg.BetaOpt; g = obj.Cfg.Gamma;
            logA = (b/g) * (1 + log(g/b));
            C = exp(logA - log(g) + gammaln(b/g));
        end
    end
    methods (Static)
        function psi = waveletFFT(w, scale, beta, gamma, side)
            ws = scale .* w;
            if strcmp(side, 'pos')
                wd = max(ws, 0); wd(w > pi) = 0;
            else
                wd = zeros(size(w)); neg = w > pi;
                wd(neg) = scale * (2*pi - w(neg));
            end
            psi = zeros(size(wd)); idx = wd > 0;
            logA = (beta/gamma) * (1 + log(gamma/beta));
            psi(idx) = exp(logA + beta .* log(wd(idx)) - wd(idx).^gamma);
        end
        function da = scaleSpacing(scales)
            scales = scales(:); Ns = length(scales);
            da = zeros(Ns, 1);
            da(1)       = abs(scales(2) - scales(1));
            da(2:end-1) = abs(scales(3:end) - scales(1:end-2)) / 2;
            da(end)     = abs(scales(end) - scales(end-1));
        end
    end
end
