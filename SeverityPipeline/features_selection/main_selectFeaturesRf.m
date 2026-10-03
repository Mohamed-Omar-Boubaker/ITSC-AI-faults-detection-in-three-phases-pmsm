projectRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(projectRoot));

cfg = SelectionConfig.rf();
FeatureSelectionPipeline(cfg, ForestRelevance(cfg.ForestTrees, cfg.ForestDepth), ...
    ForestSubsetEvaluator(cfg.NSplits, cfg.Seed, cfg.ForestTrees, cfg.ForestDepth), ...
    ForestRandomSearch(cfg.SearchIter, cfg.NSplits, cfg.Seed)).run();