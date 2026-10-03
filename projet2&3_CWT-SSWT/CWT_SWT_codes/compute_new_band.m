function [f1, f2] = compute_new_band(f0, min_scales, V)
    B_oct = min_scales / V;
    f1 = f0 * 2^(-B_oct/2);
    f2 = f0 * 2^(B_oct/2);
end
