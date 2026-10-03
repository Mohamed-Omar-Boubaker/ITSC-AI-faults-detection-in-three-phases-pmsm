classdef BandTracker
    % Adaptive dual-threshold band tracking around the target instantaneous frequency
    properties (SetAccess = immutable)
        SearchRadiusHz
        MaxExpansionHz
    end
    properties
        ShrinkageFactor = 0.125
        MaxIters = 100
    end
    methods
        function obj = BandTracker(searchRadiusHz, maxExpansionHz)
            obj.SearchRadiusHz = searchRadiusHz; obj.MaxExpansionHz = maxExpansionHz;
        end

        function [bandIdx, lowerHist] = track(obj, Tx, fAxis, IF)
            [NumFreq, NumTime] = size(Tx);
            bandIdx = cell(NumTime, 1); lowerHist = zeros(NumTime, 1);
            for idx = 1:NumTime
                curIF = IF(idx);
                [~, nominal] = min(abs(fAxis - curIF));
                searchIdx = find(fAxis >= (curIF - obj.SearchRadiusHz) & fAxis <= (curIF + obj.SearchRadiusHz));
                if isempty(searchIdx), searchIdx = nominal; end
                slice = Tx(searchIdx, idx);
                [peakAmp, off] = max(slice);
                centerIdx = searchIdx(off);
                if peakAmp < eps, bandIdx{idx} = centerIdx; continue; end

                upTh = 0.50 * peakAmp; loTh = 0.01 * peakAmp;
                for it = 1:obj.MaxIters
                    Nu = BandTracker.countPeaks(slice, upTh);
                    Nl = BandTracker.countPeaks(slice, loTh);
                    if Nu == Nl || (upTh - loTh) < 1e-7 * peakAmp, break; end
                    D = upTh - loTh;
                    upTh = upTh - obj.ShrinkageFactor * D;
                    loTh = loTh + obj.ShrinkageFactor * D;
                end
                lowerHist(idx) = loTh;

                limUp   = find(fAxis <= (curIF + obj.MaxExpansionHz), 1, 'last');
                limDown = find(fAxis >= (curIF - obj.MaxExpansionHz), 1, 'first');
                if isempty(limUp),   limUp   = NumFreq; end
                if isempty(limDown), limDown = 1; end
                up = centerIdx;
                while up < limUp && Tx(up + 1, idx) >= loTh, up = up + 1; end
                dn = centerIdx;
                while dn > limDown && Tx(dn - 1, idx) >= loTh, dn = dn - 1; end
                bandIdx{idx} = max(1, dn) : min(NumFreq, up);
            end
        end

        % Instantaneous energy of the signal reconstructed from the tracked band
        function e = reconstructEnergy(~, TxC, Tx, bandIdx, lowerHist, Cpsi)
            NumTime = length(bandIdx); e = zeros(NumTime, 1);
            for idx = 1:NumTime
                r = bandIdx{idx};
                c = TxC(r, idx);
                c(Tx(r, idx) < lowerHist(idx)) = 0;
                e(idx) = abs(sum(c, 1) / Cpsi)^2;
            end
        end
    end
    methods (Static)
        function n = countPeaks(s, th)
            above = find(s > th);
            if isempty(above), n = 0; return; end
            len = length(s);
            if len == 1, n = 1; return; end
            isPeak = false(size(s));
            for i = 2:(len-1)
                if s(i) > th && s(i) >= s(i-1) && s(i) >= s(i+1), isPeak(i) = true; end
            end
            if s(1)   > th && s(1)   > s(2),     isPeak(1)   = true; end
            if s(end) > th && s(end) > s(end-1), isPeak(end) = true; end
            n = sum(isPeak);
            if n == 0 && ~isempty(above), n = 1; end
        end
    end
end
