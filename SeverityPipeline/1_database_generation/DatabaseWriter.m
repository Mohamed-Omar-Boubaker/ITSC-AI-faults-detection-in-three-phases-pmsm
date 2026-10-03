classdef (Abstract) DatabaseWriter
    % Where the feature table is stored
    methods (Abstract)
        write(obj, featTable, names, degree, percent)
        [featTable, names, degree, percent] = loadExisting(obj)
    end
end
