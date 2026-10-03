projectRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(projectRoot));

% feature_variance.txt = the lines between the triple quotes of feats_variance.py ("name value" per line)
T = FeatureVariancePlotter.readTable(fullfile(fileparts(mfilename('fullpath')), 'feature_variance.txt'));
FeatureVariancePlotter.plot(T, 'all_98_features_variance');