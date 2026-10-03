classdef (Abstract) FeatureExtractor
    % One extraction pipeline: PhaseCurrents (already at model Fs) -> FeatureSet
    methods (Abstract)
        fs = extract(obj, currents)
    end
end
