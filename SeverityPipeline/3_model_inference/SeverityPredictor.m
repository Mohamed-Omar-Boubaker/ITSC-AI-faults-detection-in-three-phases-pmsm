classdef SeverityPredictor
    % Reproduces the training chain on one segment:
    % features -> keep trained columns -> family collapse -> variance filter -> scaling -> SVR
    properties (SetAccess = immutable)
        Bundle
        Extractors
    end
    methods
        function obj = SeverityPredictor(bundle, extractors)
            obj.Bundle = bundle; obj.Extractors = extractors;
        end

        function [sev, names] = predictSegment(obj, cur)
            b = obj.Bundle;
            sets = cellfun(@(e) e.extract(cur), obj.Extractors, 'UniformOutput', false);
            fs = FeatureSet.merge(sets);
            keep = ismember(fs.Names, b.RawFeatureNames);          % only features used in training
            if ~any(keep)
                error('SeverityPredictor:names', 'No extracted feature matches the trained feature names.');
            end
            names = fs.Names(keep);
            X  = b.Collapser.transform(fs.Values(keep), names);
            Xs = b.Scaler.transform(b.VarianceFilter.transform(X));
            sev = SvrTrainer.estimate(b.Model, Xs, true);
        end

        function res = classify(obj, segments, tStart, segmentSec)
            n = numel(segments); sev = zeros(1, n);
            for k = 1:n
                fprintf('Segment %d / %d ...\n', k, n);
                [sev(k), names] = obj.predictSegment(segments{k});
                if k == 1                                           % consistency check vs training
                    [nDiff, ex] = obj.Bundle.Collapser.countMismatch(names);
                    if nDiff > 0
                        warning('SeverityPredictor:windows', ['%d feature families have a different number of windows ' ...
                            'than in training (e.g. %s). Use a SegmentSec equal to the training sample length ' ...
                            'for features consistent with training.'], nDiff, ex);
                    end
                end
            end
            res = ClassificationResult();
            res.TStart = tStart; res.TEnd = tStart + segmentSec; res.Severity = sev;
            res.Edges = obj.Bundle.Edges; res.ClassNames = obj.Bundle.ClassNames;
            res.ClassIdx = discretize(sev, res.Edges);
        end
    end
end
