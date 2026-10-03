classdef SignalSegmenter
    % Removes the first DiscardSec seconds, then cuts the rest into SegmentSec segments
    properties (SetAccess = immutable)
        DiscardSec
        SegmentSec
    end
    methods
        function obj = SignalSegmenter(discardSec, segmentSec)
            obj.DiscardSec = discardSec; obj.SegmentSec = segmentSec;
        end
        function [segments, tStart] = split(obj, cur)
            fs = cur.Fs;
            n0 = round(obj.DiscardSec * fs);
            nS = round(obj.SegmentSec * fs);
            nTot = numel(cur.Ia);
            nSeg = floor((nTot - n0) / nS);
            if nSeg < 1
                error('SignalSegmenter:short', 'Signal too short: %.2f s left after discarding %.2f s, need >= %.2f s.', ...
                    (nTot - n0)/fs, obj.DiscardSec, obj.SegmentSec);
            end
            leftover = (nTot - n0) - nSeg * nS;
            if leftover > 0
                fprintf('Note: last %.3f s (incomplete segment) ignored.\n', leftover / fs);
            end
            segments = cell(1, nSeg); tStart = zeros(1, nSeg);
            for k = 1:nSeg
                idx = n0 + (k-1)*nS + (1:nS);
                segments{k} = PhaseCurrents(cur.Ia(idx), cur.Ib(idx), cur.Ic(idx), fs);
                tStart(k) = (idx(1) - 1) / fs;
            end
        end
    end
end
