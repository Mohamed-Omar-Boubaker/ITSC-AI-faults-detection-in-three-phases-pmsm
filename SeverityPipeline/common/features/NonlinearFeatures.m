classdef NonlinearFeatures < FeatureGroup
    % Complexity / entropy features on |x|
    methods
        function r = registry(~)
            r = { ...
                'Sample_Entropy',           @(x) NonlinearFeatures.sampleEntropy(x, 2, 0.2 * std(abs(x))); ...
                'Higuchi_Fractal_Dim',      @(x) NonlinearFeatures.higuchi(x, 10); ...
                'Lempel_Ziv_Complexity',    @NonlinearFeatures.lempelZiv; ...
                'Hilbert_Envelope_Entropy', @NonlinearFeatures.envelopeEntropy};
        end
    end
    methods (Static)
        function val = sampleEntropy(x, m, r)
            e = abs(x(:)); N = length(e);
            if nargin < 2 || isempty(m), m = 2; end
            if nargin < 3 || isempty(r), r = 0.2 * std(e, 1); end
            Xm = zeros(N - m, m); Xm1 = zeros(N - m, m + 1);
            for k = 1:m
                Xm(:, k)  = e(k : N - m + k - 1);
                Xm1(:, k) = e(k : N - m + k - 1);
            end
            Xm1(:, m+1) = e(m + 1 : N);
            B = 0; A = 0; Nm = N - m;
            for i = 1 : Nm - 1
                dm  = max(abs(Xm(i+1:end, :)  - Xm(i, :)),  [], 2);
                dm1 = max(abs(Xm1(i+1:end, :) - Xm1(i, :)), [], 2);
                B = B + sum(dm < r); A = A + sum(dm1 < r);
            end
            if A == 0 || B == 0, val = 0; else, val = -log(A / B); end
        end
        function val = higuchi(x, kmax)
            e = abs(x(:)); N = length(e);
            if nargin < 2 || isempty(kmax), kmax = 10; end
            L = zeros(kmax, 1);
            for k = 1:kmax
                Lm = zeros(k, 1);
                for m = 1:k
                    idx = m:k:N; Nm = length(idx);
                    if Nm > 1
                        Lm(m) = (sum(abs(diff(e(idx)))) * (N - 1) / (floor((N - m) / k) * k + eps)) / k;
                    end
                end
                L(k) = mean(Lm);
            end
            p = polyfit(log(1 ./ (1:kmax)'), log(L), 1); val = p(1);
        end
        function val = lempelZiv(x)
            e = abs(x(:)); N = length(e); s = e >= median(e);
            c = 1; l = 1; i = 1; k = 1; kmax = 1;
            while (i + k) <= N
                if s(i + k) == s(l + k)
                    k = k + 1;
                    if (l + k) > i
                        i = i + kmax; c = c + 1; l = 1; k = 1; kmax = 1;
                    end
                else
                    if k > kmax, kmax = k; end
                    l = l + 1;
                    if l == i
                        i = i + kmax; c = c + 1; l = 1; k = 1; kmax = 1;
                    else
                        k = 1;
                    end
                end
            end
            val = c / (N / (log2(N) + eps));
        end
        function val = envelopeEntropy(x)
            e = abs(x(:)); Et = sum(e.^2);
            if Et == 0, val = 0; return; end
            p = (e.^2) / Et; p = p(p > 0);
            val = -sum(p .* log2(p));
        end
    end
end
