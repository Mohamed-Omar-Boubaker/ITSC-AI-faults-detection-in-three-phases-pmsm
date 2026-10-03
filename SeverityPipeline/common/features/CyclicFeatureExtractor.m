classdef CyclicFeatureExtractor < FeatureExtractor
    % Pipeline 1: every WindowFeature is computed on each consecutive window
    % (cyclic analysis), then stacked feature-major (all windows of feature 1, then 2, ...)
    properties (SetAccess = immutable)
        Cfg
        WindowFeatures
    end
    methods
        function obj = CyclicFeatureExtractor(cfg, windowFeatures)
            obj.Cfg = cfg; obj.WindowFeatures = windowFeatures;
        end
        function fs = extract(obj, cur)
            c = obj.Cfg;
            t = (0:numel(cur.Ia)-1)' / c.Fs;
            [iz, ip, ineg] = SymmetricalComponents.hilbertBased(cur.Ia, cur.Ib, cur.Ic, c.Fs, c.F0);
            whole   = SignalBundle(t, cur.SpaceVector, iz, ip, ineg);
            windows = WindowSegmenter(c.WindowLength, c.Fs).segment(whole);

            nW = numel(windows); featEvo = []; base = {};
            for w = 1:nW
                sets = cellfun(@(f) f.extract(windows{w}), obj.WindowFeatures, 'UniformOutput', false);
                one  = FeatureSet.merge(sets);
                featEvo = [featEvo; one.Values];   %#ok<AGROW>
                base = one.Names;
            end

            names = cell(1, numel(featEvo));
            for k = 1:numel(base)
                names((k-1)*nW + 1 : k*nW) = FeatureSet.indexedNames(base{k}, nW);
            end
            fs = FeatureSet(names, featEvo(:)');
        end
    end
end
