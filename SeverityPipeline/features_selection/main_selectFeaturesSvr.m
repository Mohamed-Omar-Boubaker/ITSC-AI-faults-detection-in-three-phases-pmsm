projectRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(projectRoot));

cfg = SelectionConfig.svr();            % edit cfg.InFile / cfg.OutFile if needed
FeatureSelectionPipeline(cfg, StatisticalRelevance(cfg.Seed), ...
    SvrSubsetEvaluator(cfg.NSplits, cfg.Seed, cfg.BaselineC, cfg.BaselineEpsilon), ...
    SvrRandomSearch(cfg.SearchIter, cfg.NSplits, cfg.Seed)).run();