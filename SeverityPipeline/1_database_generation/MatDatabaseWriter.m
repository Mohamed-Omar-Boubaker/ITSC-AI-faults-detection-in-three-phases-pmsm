classdef MatDatabaseWriter < DatabaseWriter
    % Incremental save to a .mat file (variables: feat_table, names, degree, percent)
    properties (SetAccess = immutable)
        Path
    end
    methods
        function obj = MatDatabaseWriter(path), obj.Path = path; end
        function write(obj, featTable, names, degree, percent)
            feat_table = featTable; %#ok<NASGU>
            save(obj.Path, 'feat_table', 'names', 'degree', 'percent');
        end
        function [featTable, names, degree, percent] = loadExisting(obj)   % used to resume a run
            S = load(obj.Path);
            featTable = S.feat_table; names = S.names; degree = S.degree; percent = S.percent;
        end
    end
end
