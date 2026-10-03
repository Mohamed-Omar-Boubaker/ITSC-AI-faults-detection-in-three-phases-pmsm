classdef SampleWeighter
    % Sample weights: healthy boosted, high severity up-weighted
    properties (SetAccess = immutable)
        Thr
        HighW
        LowW
        HealthyW
    end
    methods
        function obj = SampleWeighter(cfg)
            obj.Thr = cfg.SEVERITY_THRESHOLD; obj.HighW = cfg.HIGH_WEIGHT;
            obj.LowW = cfg.LOW_WEIGHT; obj.HealthyW = cfg.HEALTHY_WEIGHT_BOOST;
        end
        function w = compute(obj, y)
            w = obj.LowW * ones(size(y));
            w(y > obj.Thr) = obj.HighW;
            w(y == 0) = obj.HealthyW;
        end
    end
end
