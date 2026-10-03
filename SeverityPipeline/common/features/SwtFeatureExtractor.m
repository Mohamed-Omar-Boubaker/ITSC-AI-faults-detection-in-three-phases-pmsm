classdef SwtFeatureExtractor < FeatureExtractor
    % Pipeline 2 (cyclic, one value per WindowSec window):
    % SWT_RecEnergy / SWT_WaveletEntropy / SWT_SubbandEnergyRatio
    properties (SetAccess = immutable)
        Cfg
        Cwt
        Tracker
    end
    methods
        function obj = SwtFeatureExtractor(cfg)
            obj.Cfg = cfg;
            obj.Cwt = MorseCwt(cfg);
            obj.Tracker = BandTracker(cfg.SearchRadiusHz, cfg.MaxExpansionHz);
        end
        function fs = extract(obj, cur)
            c = obj.Cfg; x = cur.SpaceVector; N = length(x);
            t  = (0:N-1)' / c.Fs;
            IF = c.TargetFreq * ones(N, 1);

            % 1) bilateral Morse CWT + synchrosqueezing
            [scales, f] = obj.Cwt.morseScales(N);
            da = MorseCwt.scaleSpacing(scales);
            [W, dW] = obj.Cwt.transform(x, scales);
            TxPos = SynchrosqueezeTransform.reassign(W(:,:,1), dW(:,:,1), f, scales, c.ThreshRatio, da);
            TxNeg = SynchrosqueezeTransform.reassign(W(:,:,2), dW(:,:,2), f, scales, c.ThreshRatio, da);
            clear W dW;

            Cpsi = obj.Cwt.admissibilityConstant();
            [TxAll, fAxis, TxAllC] = SynchrosqueezeTransform.buildBilateral(TxPos, TxNeg, f, Cpsi);
            clear TxPos TxNeg;

            % 2) adaptive band around the target frequency + cyclic energy metrics
            [bandIdx, lowerHist] = obj.Tracker.track(TxAll, fAxis, IF);
            WE    = SqueezeEnergyMetrics.waveletEntropy(TxAll, t, c.WindowSec);
            ratio = SqueezeEnergyMetrics.subbandEnergyRatio(TxAll, bandIdx, lowerHist, t, c.WindowSec);
            instE = obj.Tracker.reconstructEnergy(TxAllC, TxAll, bandIdx, lowerHist, Cpsi);
            rec   = SqueezeEnergyMetrics.cyclicEnergy(instE, t, c.WindowSec);

            names = [FeatureSet.indexedNames('SWT_RecEnergy', numel(rec)), ...
                     FeatureSet.indexedNames('SWT_WaveletEntropy', numel(WE)), ...
                     FeatureSet.indexedNames('SWT_SubbandEnergyRatio', numel(ratio))];
            fs = FeatureSet(names, [rec(:)', WE(:)', ratio(:)']);
        end
    end
end
