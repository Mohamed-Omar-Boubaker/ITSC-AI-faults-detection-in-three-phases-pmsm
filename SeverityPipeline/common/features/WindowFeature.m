classdef (Abstract) WindowFeature
    % A feature block computed on ONE window (SignalBundle)
    methods (Abstract)
        fs = extract(obj, window)
    end
end
