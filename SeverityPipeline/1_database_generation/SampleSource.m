classdef (Abstract) SampleSource
    % Where raw samples come from (swap HDF5 for anything else)
    methods (Abstract)
        n = count(obj)
        [currents, degree, percent] = read(obj, i)
    end
end
