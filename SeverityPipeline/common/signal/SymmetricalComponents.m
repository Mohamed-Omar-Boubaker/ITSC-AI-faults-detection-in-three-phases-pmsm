classdef SymmetricalComponents
    % Clarke space vector and Hilbert-based symmetrical components
    methods (Static)
        function is = spaceVector(ia, ib, ic)
            is = sqrt(2/3) * (ia + ib * exp(1i*2*pi/3) + ic * exp(1i*4*pi/3));
        end
        function [iZero, iPos, iNeg] = hilbertBased(ia, ib, ic, fs, f0)
            thA = angle(hilbert(ia)); thB = angle(hilbert(ib)); thC = angle(hilbert(ic));
            W = round(fs / f0);
            A = sqrt(movmean(ia.^2, W)) .* exp(1i * thA);
            B = sqrt(movmean(ib.^2, W)) .* exp(1i * thB);
            C = sqrt(movmean(ic.^2, W)) .* exp(1i * thC);
            a = exp(1i * 2 * pi / 3);
            iZero = (1/3) * (A + B + C);
            iPos  = (1/3) * (A + a   * B + a^2 * C);
            iNeg  = (1/3) * (A + a^2 * B + a   * C);
        end
    end
end
