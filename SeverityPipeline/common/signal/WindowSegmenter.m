classdef WindowSegmenter
    % Splits a SignalBundle into consecutive non-overlapping windows (cyclic analysis)
    properties (SetAccess = immutable)
        WindowLength   % seconds
        Fs
    end
    methods
        function obj = WindowSegmenter(windowLength, fs)
            obj.WindowLength = windowLength; obj.Fs = fs;
        end
        function windows = segment(obj, s)
            nWin = round(obj.WindowLength * obj.Fs);
            nW = floor(s.numSamples() / nWin);
            windows = cell(1, nW);
            for w = 1:nW
                windows{w} = s.slice((w-1)*nWin + 1 : w*nWin);
            end
        end
    end
end
