function frac = compute_energy(beta, gamma, r1, r2)
    w0 = (beta/gamma)^(1/gamma);
    k = (2*beta+1)/gamma;
    u1 = 2*(w0*r1)^gamma;
    u2 = 2*(w0*r2)^gamma;
    frac = gammainc(u2, k) - gammainc(u1, k);
end
