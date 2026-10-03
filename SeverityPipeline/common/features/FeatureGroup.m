classdef (Abstract) FeatureGroup
    % Computes labelled features from ONE signal. Subclasses only declare registry().
    methods (Abstract)
        r = registry(obj)       % n-by-2 cell: {label, @(x) value}
    end
    methods
        function fs = extract(obj, x)
            reg = obj.registry();
            vals = zeros(1, size(reg, 1));
            for k = 1:size(reg, 1), vals(k) = reg{k, 2}(x); end
            fs = FeatureSet(reg(:, 1)', vals);
        end
    end
end
