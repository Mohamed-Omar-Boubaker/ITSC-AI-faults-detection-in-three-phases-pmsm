function results = run_pipeline(sig_x, t, p)
    fprintf('\n=====================================================================\n');
    fprintf('         PRODUCTION PIPELINE REPORT\n');
    fprintf('=====================================================================\n');

    % ---- FFT detection (bilateral) ----
    [fft_freqs_all, fft_amps_all, ~, ~] = detect_frequency_FFT(sig_x, p);
    fprintf('[FFT] All peaks: %s\n', mat2str(round(fft_freqs_all,2)));

    pos_idx = fft_freqs_all > 0;
    fft_freqs_pos = fft_freqs_all(pos_idx);
    fft_amps_pos  = fft_amps_all(pos_idx);
    
    [~, max_idx] = max(fft_amps_pos); % Find index of the highest amplitude
    fft_freqs_pos = fft_freqs_pos(max_idx); % Keep only that frequency
    fft_amps_pos  = fft_amps_pos(max_idx);  % Keep only that amplitude
    
    neg_idx = fft_freqs_all < 0;
    fft_freqs_neg = fft_freqs_all(neg_idx);
    fft_amps_neg  = fft_amps_all(neg_idx);
    
    if isempty(fft_freqs_neg)
        % SAFE FALLBACK: If no negative frequency exists, treat it as 0 Hz with 0 amplitude
        fft_freqs_neg = 0; 
        fft_amps_neg  = 0;
    else
        [~, max_idx] = max(fft_amps_neg); 
        fft_freqs_neg = fft_freqs_neg(max_idx); 
        fft_amps_neg  = fft_amps_neg(max_idx);  
    end
    
    results.original_signal = sig_x;
    results.Fs = p.fs;
    results.detected_freqs = [fft_freqs_all(:), fft_amps_all(:)];

    x_rec_total = zeros(size(sig_x));

    % process only for the positive component:
    fe_k = fft_freqs_pos;
    A_k  = fft_amps_pos;
    
    % define fixed reconstruction band around targeted frequency
    fmin = max(fe_k - p.band_half_width, 0.1);
    fmax = fe_k + p.band_half_width;
    p.userBand = [fmin, fmax];
    fprintf('[BAND] Fc=%.2f Hz -> [%.2f, %.2f] Hz\n', fe_k, fmin, fmax);
    
    % validate if band conains enough scales for good reconstruction:
    if ~validate_band(fmin, fmax, 48, p.min_scales)
        [fmin, fmax] = compute_new_band(fe_k, p.min_scales, 48);
        p.userBand = [fmin, fmax];
        fprintf('[BAND] Enlarged -> [%.2f, %.2f] Hz\n', fmin, fmax);
    end
    
    % find optimal wavelet param to concentrate its energy insie the desired bandwidth
    [beta_opt, ~] = optimize_analytical_beta(fmin, fmax, p.targetEnergy, p.gamma, fe_k);
    p.beta_opt = beta_opt;
    
    % find optimal VPO : minimum overlap between wavelets - nearest scale's
    % frequency to targeted frequency fe_k - enough scales' number for
    % better reconstruction
    vpo_opt = Opt_VPO_no_overlap(beta_opt, p.gamma, fmin, fmax, fe_k);
    p.voices_per_octave = vpo_opt;
    fprintf('[VPO] Selected VPO = %d\n', vpo_opt);

    [W, f] = cwt(sig_x, "morse", p.cwt_fs, ...
                 'WaveletParameters', [p.gamma, p.gamma*beta_opt], ...
                 'VoicesPerOctave', vpo_opt);

    [x_rec_pos, band_pos] = cwt_extract_from_W(fe_k, p, W, f);

    if fft_freqs_neg == 0
        x_rec_neg = zeros(size(sig_x)); 
    else
        [x_rec_neg, ~] = cwt_extract_from_W(fft_freqs_neg, p, W, f);
    end
        
    band_neg.f_low  = -band_pos.f_high;
    band_neg.f_high = -band_pos.f_low;
    band_neg.idx    = band_pos.idx;
    band_neg.page   = 2;

    cwt_entry = struct( ...
            'freq_pos',  fft_freqs_pos, ...
            'freq_neg',  fft_freqs_neg, ...
            'amp_pos',   fft_amps_pos, ...
            'amp_neg',   fft_amps_neg, ...
            'beta',      beta_opt, ...
            'vpo',       vpo_opt, ...
            'W',         W, ...
            'f',         f, ...
            'band_pos',  band_pos, ...
            'band_neg',  band_neg, ...
            'x_rec_pos', x_rec_pos(:), ...
            'x_rec_neg', x_rec_neg(:));

    results.cwt_info = cwt_entry;

    x_rec_total = x_rec_total(:) + x_rec_pos(:) + x_rec_neg(:);

    results.x_rec_total = x_rec_total;

    fprintf('---------------------------------------------------------------------\n');
end