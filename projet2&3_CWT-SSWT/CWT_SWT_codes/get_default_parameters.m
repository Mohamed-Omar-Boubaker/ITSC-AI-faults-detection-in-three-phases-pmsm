function p = get_default_parameters()
    p.fs            = 10000;
    p.T             = 10;
    p.fe            = 100;
    p.wavelet       = 'morse';
    p.gamma         = 3;
    p.targetEnergy  = 0.99;
    p.userBand      = [];
    p.f_disp        = 1000;
    p.peak_min_dist = 5;
    p.peak_thr      = 0.01;
    p.voices_per_octave = 4;
    p.cwt_fs        = 10000;
    p.min_scales    = 2;
    p.band_half_width = 2;
end
