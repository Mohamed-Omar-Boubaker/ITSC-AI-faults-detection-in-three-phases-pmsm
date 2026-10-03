classdef Hdf5Shuffler
    % Writes a shuffled copy of the HDF5 dataset (signals and labels move together)
    properties (SetAccess = immutable)
        SrcPath
        DstPath
        Seed
    end
    methods
        function obj = Hdf5Shuffler(src, dst, seed)
            obj.SrcPath = src; obj.DstPath = dst; obj.Seed = seed;
        end
        function run(obj)
            sig = {'/current_phase_A', '/current_phase_B', '/current_phase_C'};
            sz = h5info(obj.SrcPath, sig{1}).Dataspace.Size;          % [winLen N] in MATLAB order
            winLen = sz(1); N = sz(2);
            rng(obj.Seed, 'twister'); perm = randperm(N);

            if isfile(obj.DstPath), delete(obj.DstPath); end
            for k = 1:3, h5create(obj.DstPath, sig{k}, [winLen N], 'Datatype', 'single'); end
            h5create(obj.DstPath, '/severity_percent', N, 'Datatype', 'double');
            h5create(obj.DstPath, '/severity_degree',  N, 'Datatype', 'int32');
            h5create(obj.DstPath, '/power_rating',     N, 'Datatype', 'string');

            v = h5read(obj.SrcPath, '/severity_percent'); h5write(obj.DstPath, '/severity_percent', v(perm));
            v = h5read(obj.SrcPath, '/severity_degree');  h5write(obj.DstPath, '/severity_degree',  v(perm));
            v = h5read(obj.SrcPath, '/power_rating');     h5write(obj.DstPath, '/power_rating',     v(perm));

            for new = 1:N                                              % copy signal columns in the new order
                old = perm(new);
                for k = 1:3
                    col = h5read(obj.SrcPath, sig{k}, [1 old], [winLen 1]);
                    h5write(obj.DstPath, sig{k}, col, [1 new], [winLen 1]);
                end
                if mod(new - 1, 100) == 0, fprintf('%d/%d\n', new - 1, N); end
            end
            h5writeatt(obj.DstPath, '/', 'sample_rate', h5readatt(obj.SrcPath, '/', 'sample_rate'));
            h5writeatt(obj.DstPath, '/', 'window_length_seconds', h5readatt(obj.SrcPath, '/', 'window_length_seconds'));
            fprintf('done -> %s\n', obj.DstPath);
        end
    end
end