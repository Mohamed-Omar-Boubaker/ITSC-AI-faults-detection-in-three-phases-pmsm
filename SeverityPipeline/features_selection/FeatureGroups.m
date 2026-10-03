classdef FeatureGroups
    % Columns named '<base>_wN' -> one group per base (all of its windows)
    properties (SetAccess = private)
        Bases = {}     % group names, order of first appearance
        Cols = {}      % Cols{g}: column indices ordered by window number
        Wins = {}      % Wins{g}: window numbers (same order)
        Index          % containers.Map: base -> group index
    end
    methods
        function obj = FeatureGroups(names, nCols, requiredWindows)
            tok = regexp(names, '^(.*)_w(\d+)$', 'tokens', 'once');
            bad = find(cellfun(@isempty, tok), 1);
            if ~isempty(bad), error('Column name ''%s'' does not match ''<base>_wN''.', names{bad}); end
            if numel(names) ~= nCols, error('len(names) ~= number of columns'); end
            obj.Index = containers.Map('KeyType', 'char', 'ValueType', 'double');
            for i = 1:numel(names)
                b = tok{i}{1}; w = str2double(tok{i}{2});
                if ~isKey(obj.Index, b)
                    obj.Bases{end+1} = b; obj.Index(b) = numel(obj.Bases);
                    obj.Cols{end+1} = []; obj.Wins{end+1} = [];
                end
                g = obj.Index(b);
                obj.Cols{g}(end+1) = i; obj.Wins{g}(end+1) = w;
            end
            for g = 1:numel(obj.Bases)
                [obj.Wins{g}, o] = sort(obj.Wins{g}); obj.Cols{g} = obj.Cols{g}(o);
                if ~isempty(requiredWindows) && ~isequal(obj.Wins{g}, 1:requiredWindows)
                    error('Group ''%s'' has windows [%s], expected 1-%d.', obj.Bases{g}, num2str(obj.Wins{g}), requiredWindows);
                end
            end
        end
        function n = numGroups(obj), n = numel(obj.Bases); end
        function tf = has(obj, b), tf = isKey(obj.Index, b); end
        function c = firstCol(obj, b), c = min(obj.Cols{obj.Index(b)}); end
        function c = cols(obj, bases)            % all columns of the given groups, in the given order
            c = [];
            for k = 1:numel(bases), c = [c, obj.Cols{obj.Index(bases{k})}]; end %#ok<AGROW>
        end
    end
end