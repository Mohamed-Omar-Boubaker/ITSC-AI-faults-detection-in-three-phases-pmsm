function ok = validate_band(f1, f2, V, min_scales)
    N = floor(V * log2(f2/f1));
    ok = (N >= min_scales);
    if ~ok
        fprintf('[WARN] Band too narrow (scales=%d < %d)\n', N, min_scales);
    end
end
