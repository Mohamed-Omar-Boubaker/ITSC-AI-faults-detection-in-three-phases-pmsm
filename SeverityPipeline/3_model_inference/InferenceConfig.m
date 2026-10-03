classdef InferenceConfig
    % Settings of the inference run (edit them in main_classifySignal.m)
    properties
        BundlePath = 'C:\Users\moham\OneDrive\Desktop\SeverityPipeline\svr_model_bundle.mat'
        SignalFile = 'C:\Users\moham\OneDrive\Desktop\SeverityPipeline\itsc_03pct.mat'
        VarNames = {'ia', 'ib', 'ic'}   % 3 variables (phases A,B,C) or 1 variable (N-by-3 matrix)
        FsSignal = []                   % sampling rate of the file; [] -> read variable 'fs' in the file
        FsModel = 10000                 % rate used by the feature extractors (training rate)
        DiscardSec = 2                  % initial seconds removed (transient)
        SegmentSec = 2                  % length of each classified segment
    end
end
