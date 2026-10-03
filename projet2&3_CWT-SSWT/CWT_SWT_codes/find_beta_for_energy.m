function beta = find_beta_for_energy(r1, r2, gamma, target_energy)
    beta_low = 0.1; beta_high = 1.0; max_beta = 300;
    frac_ceiling = compute_energy(max_beta, gamma, r1, r2);
    if frac_ceiling < target_energy
        target_energy = 0.97 * frac_ceiling;
    end
    f_low  = compute_energy(beta_low,  gamma, r1, r2) - target_energy;
    f_high = compute_energy(beta_high, gamma, r1, r2) - target_energy;
    while f_high <= 0 && beta_high < max_beta
        beta_high = beta_high * 2;
        f_high = compute_energy(beta_high, gamma, r1, r2) - target_energy;
    end
    if f_low < 0 && f_high > 0
        beta = fzero(@(b) compute_energy(b, gamma, r1, r2) - target_energy, [beta_low, beta_high]);
    else
        beta = fminbnd(@(b) abs(compute_energy(b, gamma, r1, r2) - target_energy), 0.1, min(max_beta,100));
    end
end
