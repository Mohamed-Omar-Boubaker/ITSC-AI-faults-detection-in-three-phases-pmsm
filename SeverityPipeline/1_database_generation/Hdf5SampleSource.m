classdef Hdf5SampleSource < SampleSource
    % Reads phase currents + severity labels from the PMSM ITSC HDF5 file
    properties (SetAccess = immutable)
        File
        WinLen
        NumSamples
        Fs
    end
    methods
        function obj = Hdf5SampleSource(file, fs)
            info = h5info(file, '/current_phase_A');
            obj.File = file; obj.Fs = fs;
            obj.WinLen = info.Dataspace.Size(1);
            obj.NumSamples = info.Dataspace.Size(2);
        end
        function n = count(obj), n = obj.NumSamples; end
        function [cur, degree, percent] = read(obj, i)
            rd = @(name) double(h5read(obj.File, name, [1 i], [obj.WinLen 1]));
            cur = PhaseCurrents(rd('/current_phase_A'), rd('/current_phase_B'), rd('/current_phase_C'), obj.Fs);
            degree  = h5read(obj.File, '/severity_degree',  i, 1);
            percent = h5read(obj.File, '/severity_percent', i, 1);
        end
    end
end
