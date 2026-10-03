classdef FeatureSet
    % Names + values of a feature vector (both row vectors)
    properties
        Names = {}
        Values = []
    end
    methods
        function obj = FeatureSet(names, values)
            if nargin > 0, obj.Names = names; obj.Values = values; end
        end
        function obj = withPrefix(obj, prefix)
            obj.Names = strcat(prefix, obj.Names);
        end
    end
    methods (Static)
        function out = merge(sets)              % concatenate a cell array of FeatureSet
            names = {}; vals = [];
            for k = 1:numel(sets)
                names = [names, sets{k}.Names];   %#ok<AGROW>
                vals  = [vals,  sets{k}.Values];  %#ok<AGROW>
            end
            out = FeatureSet(names, vals);
        end
        function n = indexedNames(base, count)  % base_w1 ... base_wN
            n = arrayfun(@(k) sprintf('%s_w%d', base, k), 1:count, 'UniformOutput', false);
        end
    end
end
