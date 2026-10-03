function [beta_opt, f0_opt] = optimize_analytical_beta(f1, f2, target_energy, gamma, f0_fixed)
    f0_opt = abs(f0_fixed);
    r1 = f1 / f0_opt;
    r2 = f2 / f0_opt;
    beta_opt = find_beta_for_energy(r1, r2, gamma, target_energy);
end
