classdef SignalBundle
    % Time axis + complex signals (space vector, zero/pos/neg sequences)
    properties
        T
        Is
        IZero
        IPos
        INeg
    end
    methods
        function obj = SignalBundle(t, is, iZero, iPos, iNeg)
            obj.T = t; obj.Is = is; obj.IZero = iZero; obj.IPos = iPos; obj.INeg = iNeg;
        end
        function n = numSamples(obj), n = numel(obj.IPos); end
        function w = slice(obj, idx)
            w = SignalBundle(obj.T(idx), obj.Is(idx), obj.IZero(idx), obj.IPos(idx), obj.INeg(idx));
        end
    end
end
