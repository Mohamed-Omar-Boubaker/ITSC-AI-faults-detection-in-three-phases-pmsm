classdef ModelBundle
    % Trained model + the fitted preprocessing chain, loaded from svr_model_bundle.mat
    properties (SetAccess = private)
        Model
        Collapser
        VarianceFilter
        Scaler
        RawFeatureNames
        Edges
        ClassNames
    end
    methods (Static)
        function b = fromFile(path)
            S = load(path);
            need = {'model', 'collapser', 'variance_filter', 'scaler', 'raw_feature_names', 'edges_final', 'class_names_final'};
            miss = need(~isfield(S, need));
            if ~isempty(miss)
                error('ModelBundle:old', ['Bundle lacks: %s. It was saved by an older script - ' ...
                    'retrain with 2_svr_training/main_trainSvr.m.'], strjoin(miss, ', '));
            end
            b = ModelBundle();
            b.Model = S.model; b.Collapser = S.collapser; b.VarianceFilter = S.variance_filter;
            b.Scaler = S.scaler; b.RawFeatureNames = S.raw_feature_names;
            b.Edges = S.edges_final; b.ClassNames = S.class_names_final;
        end
    end
end
