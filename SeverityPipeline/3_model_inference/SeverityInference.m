classdef SeverityInference
    % Orchestrator: load model + signal -> resample -> drop first seconds -> segment -> classify -> plot
    properties (SetAccess = immutable)
        Cfg
    end
    methods
        function obj = SeverityInference(cfg), obj.Cfg = cfg; end

        function res = run(obj)
            c = obj.Cfg;
            bundle     = ModelBundle.fromFile(c.BundlePath);
            extractors = FeatureExtractorFactory.createDefault(c.FsModel);   % same features as training
            predictor  = SeverityPredictor(bundle, extractors);

            cur = MatSignalLoader.read(c).resampleTo(c.FsModel);             % same rate as training
            fprintf('Signal: %.2f s at %g Hz (after resampling)\n', numel(cur.Ia)/c.FsModel, c.FsModel);

            [segments, tStart] = SignalSegmenter(c.DiscardSec, c.SegmentSec).split(cur);
            fprintf('%d segment(s) of %.1f s after discarding the first %.1f s\n', numel(segments), c.SegmentSec, c.DiscardSec);

            res = predictor.classify(segments, tStart, c.SegmentSec);
            res.print();
            SeverityEvolutionPlotter.plot(res);
        end
    end
end
