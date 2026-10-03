classdef SvrConfig
    % Everything tunable (same names / values as the original CFG struct)
    properties
        % -- severity / weighting --
        SEVERITY_THRESHOLD = 6
        HEALTHY_WEIGHT_BOOST = 20
        HIGH_WEIGHT = 2.5
        LOW_WEIGHT = 1.0
        HEALTHY_THRESHOLD = 1.3
        % -- targets / SVR search --
        TARGET_SV_RATIO = 0.3
        TARGET_ACCURACY = 0.9
        SV_PENALTY_LAMBDA = 60
        BOX_RANGE = [0.11, 5000]
        KS_RANGE = [0.05, 500]
        EPS_RANGE = [0.002, 3]
        VERBOSE_TRIALS = true
        % -- nested CV (optional) --
        RUN_NESTED_CV = false
        NESTED_OUTER_FOLDS = 5
        NESTED_OUTER_REPEATS = 3
        NESTED_INNER_FOLDS = 5
        NESTED_INNER_REPEATS = 1
        NESTED_INNER_MAX_EVALS = 25
        % -- final search --
        FINAL_SEARCH_FOLDS = 5
        FINAL_SEARCH_REPEATS = 3
        FINAL_SEARCH_MAX_EVALS = 60
        % -- multiclass boundary search (DP) --
        MIN_TOTAL_CLASSES = 4
        MAX_TOTAL_CLASSES = 6
        MIN_CLASS_FRAC = 0.05
        MAX_CLASS_FRAC = 0.75
        COMPLEXITY_LAMBDA = 0.08
        BOUNDARY_CV_FOLDS = 5
        BOUNDARY_CV_REPEATS = 4
        SCORE_WEIGHTS = [0.60, 0.25, 0.15]   % [accuracy, macroF1, macroMCC]
        USE_BALANCED_ACC = false
        N_BOOT = 500
        N_CANDIDATE_GRID = 300
        TARGET_CLASS_WIDTH = 6
        WIDTH_BONUS_WEIGHT = 0.2
        % -- files --
        DATA_PATH = 'C:\Users\moham\OneDrive\Desktop\ITSC AI faults detection\SeverityPipeline\features_selection\filtered_DB_SVR_selected.mat'
        SAVE_PATH = 'C:\Users\moham\OneDrive\Desktop\SeverityPipeline\svr_model_bundle.mat'
    end
    methods
        function s = toStruct(obj)      % plain struct saved in the bundle (like the old CFG)
            for p = properties(obj)'
                s.(p{1}) = obj.(p{1});
            end
        end
    end
end
