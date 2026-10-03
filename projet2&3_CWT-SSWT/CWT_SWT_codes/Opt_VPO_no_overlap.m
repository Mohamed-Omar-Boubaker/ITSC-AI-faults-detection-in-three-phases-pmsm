function vpo_opt = Opt_VPO_no_overlap(beta, gamma, fmin, fmax, fe)
    BW = abs(fe / sqrt(2*beta*gamma));
    raw_v_min = 1 / log2(fmax/fmin);
    V_min = ceil(raw_v_min);
    if mod(V_min,2) ~= 0, V_min = V_min+1; end
    V_min = max(V_min, 4);

    for V = V_min:2:48
        N = floor(V * log2(fmax/fmin) + 1);
        f = fmax * 2.^(-(0:N-1)/V);
        [~, idx] = min(abs(f - abs(fe)));
        if idx == 1 || idx == N, continue; end
        if (abs(f(idx)-f(idx-1)) >= BW) && (abs(f(idx+1)-f(idx)) >= BW)
            vpo_opt = V;
            return;
        end
    end
    vpo_opt = V_min;  % fallback
end
