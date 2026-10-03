classdef TrialLog < handle
    % Log of every hyper-parameter trial (replaces the old global TRIAL_LOG)
    properties
        Entries = struct('BoxConstraint', {}, 'KernelScale', {}, 'Epsilon', {}, ...
                         'rmse', {}, 'sv_ratio', {}, 'cost', {})
    end
    methods
        function add(obj, e), obj.Entries(end+1) = e; end
    end
end
