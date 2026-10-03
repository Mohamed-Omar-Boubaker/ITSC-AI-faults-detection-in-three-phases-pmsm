classdef FeatureExtractorFactory
    % Single place that defines WHICH features are extracted.
    % Used by the database builder AND the inference script -> identical features.
    methods (Static)
        function extractors = createDefault(fsTarget)
            fcfg = FeatureConfig(); fcfg.Fs = fsTarget;
            scfg = SwtConfig();     scfg.Fs = fsTarget;

            time = TimeFeatures();
            freq = FrequencyFeatures(fcfg.Fs, fcfg.F0, fcfg.H);
            geom = GeometricFeatures();
            nonl = NonlinearFeatures();

            % Per-window (cyclic) feature blocks, in the original order
            windowFeatures = { ...
                ChannelFeature('I1_',  'IPos', time), ChannelFeature('I1_',  'IPos', freq), ...
                ChannelFeature('I2_',  'INeg', time), ChannelFeature('I2_',  'INeg', freq), ...
                UnbalanceFactorFeature(), ...
                ChannelFeature('CSV_', 'Is',   time), ChannelFeature('CSV_', 'Is',   freq), ...
                ChannelFeature('CSV_', 'Is',   geom), ChannelFeature('CSV_', 'Is',   nonl)};

            extractors = {CyclicFeatureExtractor(fcfg, windowFeatures), SwtFeatureExtractor(scfg)};
        end
    end
end
