classdef MatSignalLoader
    % Reads the three phase currents (and sampling rate) from a .mat file
    methods (Static)
        function cur = read(cfg)
            S = load(cfg.SignalFile);
            names = cfg.VarNames;
            missing = names(~isfield(S, names));
            if ~isempty(missing)
                error('MatSignalLoader:var', 'Variable(s) not found in %s: %s (found: %s)', ...
                    cfg.SignalFile, strjoin(missing, ', '), strjoin(fieldnames(S)', ', '));
            end
            if numel(names) == 3
                ia = S.(names{1}); ib = S.(names{2}); ic = S.(names{3});
            else
                M = S.(names{1});
                if size(M, 1) == 3 && size(M, 2) ~= 3, M = M.'; end
                ia = M(:, 1); ib = M(:, 2); ic = M(:, 3);
            end
            fs = cfg.FsSignal;
            if isempty(fs)
                if isfield(S, 'fs'), fs = double(S.fs);
                else, error('MatSignalLoader:fs', 'Set cfg.FsSignal (no variable ''fs'' in the file).');
                end
            end
            cur = PhaseCurrents(double(ia), double(ib), double(ic), fs);
        end
    end
end
