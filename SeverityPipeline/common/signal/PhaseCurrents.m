classdef PhaseCurrents
    % Three phase currents (column vectors) sampled at Fs
    properties (SetAccess = immutable)
        Ia
        Ib
        Ic
        Fs
    end
    properties (Dependent)
        SpaceVector   % complex Clarke space vector
    end
    methods
        function obj = PhaseCurrents(ia, ib, ic, fs)
            obj.Ia = ia(:); obj.Ib = ib(:); obj.Ic = ic(:); obj.Fs = fs;
        end
        function v = get.SpaceVector(obj)
            v = SymmetricalComponents.spaceVector(obj.Ia, obj.Ib, obj.Ic);
        end
        function out = resampleTo(obj, fsTarget)
            out = PhaseCurrents(SignalResampler.apply(obj.Ia, obj.Fs, fsTarget), ...
                                SignalResampler.apply(obj.Ib, obj.Fs, fsTarget), ...
                                SignalResampler.apply(obj.Ic, obj.Fs, fsTarget), fsTarget);
        end
    end
end
