%% ========================================================================
%  RUN THIS FILE TO CLASSIFY A SIGNAL WITH THE TRAINED MODEL
%  .mat signal -> drop first 2 s -> 2 s segments -> severity + class per segment -> plots
%% ========================================================================
projectRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(projectRoot));

cfg = InferenceConfig();
cfg.BundlePath = 'C:\Users\moham\OneDrive\Desktop\SeverityPipeline\svr_model_bundle.mat';
cfg.SignalFile = 'C:\Users\moham\OneDrive\Desktop\SeverityPipeline\test_30pct.mat';        % .mat file containing the phase currents
cfg.VarNames = {'ia','ib','ic'};     % or one N-by-3 variable, e.g. {'i_abc'}
cfg.FsSignal   = [];                     % [] -> use variable 'fs' stored in the file
cfg.DiscardSec = 2;                      % transient removed at the start
cfg.SegmentSec = 2;                      % segment length

result = SeverityInference(cfg).run();   % result.Severity, result.ClassIdx, result.ClassNames
