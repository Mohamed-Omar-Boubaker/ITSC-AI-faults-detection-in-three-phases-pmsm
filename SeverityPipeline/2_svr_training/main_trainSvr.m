%% ========================================================================
%  RUN THIS FILE TO TRAIN THE SVR SEVERITY MODEL AND SAVE THE MODEL BUNDLE
%  (input: filtered feature database  ->  output: svr_model_bundle.mat)
%% ========================================================================
projectRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(projectRoot));

cfg = SvrConfig();      % edit paths / parameters here, e.g.
% cfg.DATA_PATH = 'C:\...\filtered_DB_SVR_selected.mat';
% cfg.SAVE_PATH = 'C:\...\svr_model_bundle.mat';
% cfg.RUN_NESTED_CV = true;

SeverityExperiment(cfg).run();
