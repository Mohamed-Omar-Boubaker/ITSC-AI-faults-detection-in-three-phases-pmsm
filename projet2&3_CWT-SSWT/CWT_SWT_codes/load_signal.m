function [sig_x, t, p] = load_signal(filepath, p)
    if exist(filepath, 'file') == 2
        loaded = load(filepath);
        sig_x = loaded.signals;
        t = loaded.t;
        fprintf('data loaded.\n');
    else
        fprintf('[WARNING] File not found. Generating dummy signal.\n');
    end
    p.N = numel(t);
end