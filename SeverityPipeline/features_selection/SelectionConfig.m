classdef SelectionConfig
    % All selection settings. Presets: SelectionConfig.svr() / SelectionConfig.rf()
    properties
        InFile  = 'c:/Users/moham/OneDrive/Desktop/production/Full_DB_ordered.mat'
        OutFile = 'C:\Users\moham\OneDrive\Desktop\production\filtered_DB_SVR_selected.mat'
        Seed = 42
        ProtectedGroups = {'SWT_RecEnergy', 'SWT_WaveletEntropy', 'SWT_SubbandEnergyRatio', 'CSV_Eccentricity'}
        RedundancyThreshold = 0.95
        RmseTolerance = 0.05
        RelevancePercentile = 60
        MinRetainedGroups = 20
        MinFinalGroups = 30
        NSplits = 5
        Patience = 10
        CandidatesPerIter = Inf            % Inf = test every remaining group
        RequiredWindows = []               % [] = any number of windows per group
        RecomputeRelevancePerFold = true   % stability analysis redoes relevance + redundancy per fold
        SaveCleaned = true                 % true: cleaned X,y (Inf->NaN, NaN targets dropped); false: raw table + degree
        BaselineC = 496
        BaselineEpsilon = 0.9
        SearchIter = 40
        ForestTrees = 60
        ForestDepth = 10
    end
    methods (Static)
        function c = svr(), c = SelectionConfig(); end
        function c = rf()
            c = SelectionConfig();
            c.OutFile = 'C:\Users\moham\OneDrive\Desktop\production\filtered_DB_RF_selected.mat';
            c.RmseTolerance = 0.08; c.MinRetainedGroups = 30; c.NSplits = 3; c.Patience = 8;
            c.CandidatesPerIter = 10; c.RequiredWindows = 10; c.RecomputeRelevancePerFold = false;
            c.SaveCleaned = false; c.SearchIter = 15;
        end
    end
end