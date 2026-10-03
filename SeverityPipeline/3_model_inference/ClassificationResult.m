classdef ClassificationResult
    % Severity + class of every segment
    properties
        TStart       % segment start time in the original signal (s)
        TEnd
        Severity     % predicted severity (%)
        ClassIdx     % class index (1 = healthy)
        ClassNames
        Edges
    end
    methods
        function print(obj)
            fprintf('\n%-4s %-16s %-14s %s\n', 'Seg', 'Time (s)', 'Severity (%)', 'Class');
            for k = 1:numel(obj.Severity)
                fprintf('%-4d %6.2f - %-7.2f %-14.3f %s\n', k, obj.TStart(k), obj.TEnd(k), ...
                    obj.Severity(k), obj.ClassNames{obj.ClassIdx(k)});
            end
        end
    end
end
