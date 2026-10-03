classdef SignalResampler
    % Rational-factor resampling (same call as the original pipeline)
    methods (Static)
        function y = apply(x, fsOrig, fsTarget)
            [p, q] = rat(fsTarget / fsOrig);
            y = resample(x(:), p, q);
        end
    end
end
