function SSTFT()

    config = get_default_parameters();

    % Define the signal sources to process sequentially
    signal_sources = {
        struct('name', 'Loaded Signal',    'loader', @(cfg) load_signal(cfg)), ...
        struct('name', 'Synthetic Signal', 'loader', @(cfg) synthetic_signal()) ...
    };

    for s = 1:length(signal_sources)
        sig_name   = signal_sources{s}.name;
        sig_loader = signal_sources{s}.loader;

        fprintf('\n========================================================\n');
        fprintf('  PROCESSING SIGNAL %d/%d: %s\n', s, length(signal_sources), sig_name);
        fprintf('========================================================\n');

        % Load/Generate signal
        [t, x] = sig_loader(config);

        % FREQUENCY DETECTION
        [f_pos, A_pos, f_neg, A_neg, fax, fft_orig] = detect_frequency_FFT(x, config);

        fprintf('\n[%s] Detected frequency components:\n', sig_name);
        fprintf('Positive frequency: %.4f Hz | Amplitude: %.4f\n', f_pos, A_pos);
        fprintf('Negative frequency: %.4f Hz | Amplitude: %.4f\n', f_neg, A_neg);

        if config.plot_signal
            plot_signal(t, x, config);
        end
        if config.plot_fft
            plot_fft(fax, fft_orig, config, f_pos, f_neg);
        end

        % METHOD 1: Single-Window Entropy Optimization
        [win_opt, config_opt, ~] = optimize_gaussian_window(x, config);
        config_opt.C_g = sum(win_opt);
        config_opt.g_0 = win_opt(config_opt.Lh + 1);

        [STFT_single, dstft_single, f_stft, t_stft] = compute_stft( ...
            x, win_opt, config_opt.gaussian_derivative, config.fs, config.nfft, config.hop);

        STFT_single_phys = abs(STFT_single) / config_opt.C_g;

        SSTFT_single_raw  = compute_sst(STFT_single, dstft_single, f_stft, 1e-4);
        SSTFT_single_phys = single(abs(SSTFT_single_raw) / (config_opt.g_0 * config.nfft));

        % METHOD 2: Multi-Window Gauss Fusion (Memory-Optimized)
        [STFT_multi_phys, multi_info, STFT_multi_coh, dstft_multi_coh, C_g_multi, g0_multi] = ...
            compute_multi_window_stft(x, config);

        SSTFT_multi      = compute_sst(STFT_multi_coh, dstft_multi_coh, f_stft, 1e-4);
        SSTFT_multi_phys = single(abs(SSTFT_multi) / (g0_multi * config.nfft));

        % PERFORMANCE & QUANTITATIVE COMPARISON
        entropy_single = calculate_renyi_entropy(STFT_single_phys, config.renyi_alpha);
        entropy_multi  = calculate_renyi_entropy(STFT_multi_phys, config.renyi_alpha);

        zero_mask_single = STFT_single_phys < (1e-3 * max(STFT_single_phys(:)));
        zero_mask_multi  = STFT_multi_phys < (1e-3 * max(STFT_multi_phys(:)));
        ratio_zero_single = mean(zero_mask_single(:)) * 100;
        ratio_zero_multi  = mean(zero_mask_multi(:)) * 100;

        fprintf('\n======================================================\n');
        fprintf('  METHOD COMPARISON SUMMARY (%s)\n', sig_name);
        fprintf('========================================================\n');
        fprintf('Metric                     | Single-Window | Multi-Window\n');
        fprintf('----------------------------+---------------+----------------\n');
        fprintf('Rényi Entropy (Alpha=%.1f)  |    %13.4f     |    %10.4f\n', config.renyi_alpha, entropy_single, entropy_multi);
        fprintf('Low Coverage / Dead Zone %% |    %12.2f%%   |    %9.2f%%\n', ratio_zero_single, ratio_zero_multi);
        fprintf('Number of Windows Used     | %13d       |    %10d\n', 1, multi_info.N_angles);
        fprintf('=============================================================\n\n');

        % PLOTTING: spectrograms
        if config.plot_stft
            plot_stft(t_stft, f_stft, STFT_single_phys, config, f_pos, f_neg);
            title(sprintf('[%s] Single-Window (Entropy Optimized) STFT Amplitude', sig_name));

            plot_stft(t_stft, f_stft, STFT_multi_phys, config, f_pos, f_neg);
            title(sprintf('[%s] Multi-Window Gauss Fusion STFT Amplitude (N = %d angles)', sig_name, config.multi_N_angles));

            plot_stft(t_stft, f_stft, SSTFT_single_phys, config, f_pos, f_neg);
            title(sprintf('[%s] SST single window method', sig_name));

            plot_stft(t_stft, f_stft, SSTFT_multi_phys, config, f_pos, f_neg);
            title(sprintf('[%s] SST multi window method', sig_name));
        end

        % SIGNAL RECONSTRUCTION 
        sample_indices = 1:config.hop:length(x);
        x_ref = x(sample_indices);

        x_rec_stft_single  = sum(STFT_single, 1).'      / (config_opt.g_0 * config.nfft);
        x_rec_sstft_single = sum(SSTFT_single_raw, 1).' / (config_opt.g_0 * config.nfft);

        x_rec_stft_multi   = sum(STFT_multi_coh, 1).'   / (g0_multi * config.nfft);
        x_rec_sstft_multi  = sum(SSTFT_multi, 1).'      / (g0_multi * config.nfft);

        max_val = max(abs([x_ref(:); x_rec_stft_single(:); x_rec_sstft_single(:); ...
                           x_rec_stft_multi(:); x_rec_sstft_multi(:)]));
        ylims = [-max_val, max_val] * 1.05;

        plot_reconstruction(t_stft, x_ref, x_rec_stft_single,  sprintf('[%s] Single-Window: STFT Reconstruction', sig_name), ylims);
        plot_reconstruction(t_stft, x_ref, x_rec_sstft_single, sprintf('[%s] Single-Window: SST Reconstruction', sig_name),  ylims);
        plot_reconstruction(t_stft, x_ref, x_rec_stft_multi,   sprintf('[%s] Multi-Window: STFT Reconstruction', sig_name),  ylims);
        plot_reconstruction(t_stft, x_ref, x_rec_sstft_multi,  sprintf('[%s] Multi-Window: SST Reconstruction', sig_name),   ylims);

        % RECONSTRUCTION METRICS 
        Lh_single = config_opt.Lh;
        Lh_multi  = floor((config.nfft - 1) / 2);

        valid_single = (sample_indices > Lh_single) & (sample_indices <= (length(x) - Lh_single));
        valid_multi  = (sample_indices > Lh_multi)  & (sample_indices <= (length(x) - Lh_multi));

        m_stft_single  = compute_recon_metrics(x_ref(valid_single), x_rec_stft_single(valid_single));
        m_sstft_single = compute_recon_metrics(x_ref(valid_single), x_rec_sstft_single(valid_single));
        m_stft_multi   = compute_recon_metrics(x_ref(valid_multi),  x_rec_stft_multi(valid_multi));
        m_sstft_multi  = compute_recon_metrics(x_ref(valid_multi),  x_rec_sstft_multi(valid_multi));

        fprintf('\n===================================================================\n');
        fprintf('   RECONSTRUCTION METRICS SUMMARY (%s)\n', sig_name);
        fprintf('===================================================================\n');
        fprintf('Metric          | Single-STFT | Single-SST |  Multi-STFT | Multi-SST\n');
        fprintf('----------------+-------------+------------+-------------+----------\n');
        fprintf('RMSE            | %11.4e | %10.4e | %11.4e | %8.4e\n', ...
            m_stft_single.RMSE, m_sstft_single.RMSE, m_stft_multi.RMSE, m_sstft_multi.RMSE);
        fprintf('NRMSE           | %11.4f | %10.4f | %11.4f | %8.4f\n', ...
            m_stft_single.NRMSE, m_sstft_single.NRMSE, m_stft_multi.NRMSE, m_sstft_multi.NRMSE);
        fprintf('Relative err %% | %10.2f%% | %9.2f%% | %10.2f%% | %7.2f%%\n', ...
            m_stft_single.RelErr_pct, m_sstft_single.RelErr_pct, m_stft_multi.RelErr_pct, m_sstft_multi.RelErr_pct);
        fprintf('Correlation     | %11.4f | %10.4f | %11.4f | %8.4f\n', ...
            m_stft_single.Corr, m_sstft_single.Corr, m_stft_multi.Corr, m_sstft_multi.Corr);
        fprintf('===================================================================\n\n');

        % ENERGY ANALYSIS (positive/negative frequencies)
        fprintf('\n========================================================\n');
        fprintf('   ENERGY ANALYSIS (%s)\n', sig_name);
        fprintf('========================================================\n');

        E_orig = sum(abs(x_ref).^2);
        E_rec_stft_single  = sum(abs(x_rec_stft_single).^2);
        E_rec_sstft_single = sum(abs(x_rec_sstft_single).^2);
        E_rec_stft_multi   = sum(abs(x_rec_stft_multi).^2);
        E_rec_sstft_multi  = sum(abs(x_rec_sstft_multi).^2);

        fprintf('\nTime-domain energy (sum |x|^2):\n');
        fprintf('  Original signal (sampled)    : %.6e\n', E_orig);
        fprintf('  Reconstructed (STFT single)  : %.6e\n', E_rec_stft_single);
        fprintf('  Reconstructed (SST single)   : %.6e\n', E_rec_sstft_single);
        fprintf('  Reconstructed (STFT multi)   : %.6e\n', E_rec_stft_multi);
        fprintf('  Reconstructed (SST multi)    : %.6e\n', E_rec_sstft_multi);

        f_pos_idx = f_stft > 0;
        f_neg_idx = f_stft < 0;

        E_STFT_single_pos  = sum(abs(STFT_single(f_pos_idx, :)).^2, 'all');
        E_STFT_single_neg  = sum(abs(STFT_single(f_neg_idx, :)).^2, 'all');
        E_SSTFT_single_pos = sum(abs(SSTFT_single_raw(f_pos_idx, :)).^2, 'all');
        E_SSTFT_single_neg = sum(abs(SSTFT_single_raw(f_neg_idx, :)).^2, 'all');

        E_STFT_multi_pos   = sum(abs(STFT_multi_coh(f_pos_idx, :)).^2, 'all');
        E_STFT_multi_neg   = sum(abs(STFT_multi_coh(f_neg_idx, :)).^2, 'all');
        E_SSTFT_multi_pos  = sum(abs(SSTFT_multi(f_pos_idx, :)).^2, 'all');
        E_SSTFT_multi_neg  = sum(abs(SSTFT_multi(f_neg_idx, :)).^2, 'all');

        fprintf('\nFrequency-domain energy (positive/negative bands):\n');
        fprintf('Method        | Transform | Positive Energy | Negative Energy\n');
        fprintf('--------------+-----------+-----------------+---------------\n');
        fprintf('Single-Window | STFT      | %15.6e | %15.6e\n', E_STFT_single_pos, E_STFT_single_neg);
        fprintf('Single-Window | SSTFT     | %15.6e | %15.6e\n', E_SSTFT_single_pos, E_SSTFT_single_neg);
        fprintf('Multi-Window  | STFT      | %15.6e | %15.6e\n', E_STFT_multi_pos, E_STFT_multi_neg);
        fprintf('Multi-Window  | SSTFT     | %15.6e | %15.6e\n', E_SSTFT_multi_pos, E_SSTFT_multi_neg);
        fprintf('========================================================\n\n');

        % Clear memory at the end of each iteration
        clear STFT_single dstft_single SSTFT_single_raw STFT_multi_coh dstft_multi_coh SSTFT_multi;
    end

end

%% MULTI-WINDOW GAUSS FUSION STFT (FrFT-CORRECT, SCALE-DECOUPLED, MEMORY EFFICIENT)
function [STFT_fused, info, STFT_fused_coherent, dstft_fused_coherent, C_g_eff, g0_eff] = compute_multi_window_stft(x, config)

    N_angles  = config.multi_N_angles;
    a_eff     = config.multi_shape_a;     % shape: controls diversity between orientations
    sigma_max = config.multi_sigma_max;   % scale: physical width of the widest orientation
    fs   = config.fs;
    nfft = config.nfft;
    hop  = config.hop;

    % Time-rescale so the widest orientation lands exactly on sigma_max,
    % independent of a_eff. This is what lets us keep real angle-to-angle
    % diversity (a_eff != 1) while staying inside the nfft budget.
    sigma_norm_0   = 1 / sqrt(2*pi*a_eff);
    sigma_norm_pi2 = sqrt(a_eff / (2*pi));
    T_ref = sigma_max / max(sigma_norm_0, sigma_norm_pi2);

    angles = (1:N_angles) * (pi / N_angles);
    Lh_max = floor((nfft - 1) / 2);

    N_time = length(1:hop:length(x));
    STFT_fused_sum = zeros(nfft, N_time, 'single');
    n_included = 0;

    % Common length = 2*Lh_max+1 so every g_rot (length <= that, by
    % construction, since Lh = min(Lh_full, Lh_max)) fits centered inside it.
    g_eff  = zeros(2*Lh_max + 1, 1);
    dg_eff = zeros(2*Lh_max + 1, 1);

    fprintf('\nMulti-window Gauss fusion (FrFT-correct, scale-decoupled)\n');

    for j = 1:N_angles
        alpha = angles(j);

        if abs(sin(alpha)) < 1e-10
            b_j = a_eff;
            A_j = a_eff^(1/4);
        else
            cotA = cot(alpha);
            p    = a_eff - 1i*cotA;
            b_j  = 1/(p*sin(alpha)^2) - 1i*cotA;
            c_a  = sqrt(1 - 1i*cotA);
            A_j  = c_a * a_eff^(1/4) * p^(-1/2);
        end

        b_j = b_j / T_ref^2;   % apply the physical scale

        sigma_eff = 1 / sqrt(2*pi*real(b_j));
        Lh_full   = ceil(4*sigma_eff*fs);
        Lh        = max(min(Lh_full, Lh_max), 1);

        fprintf('  [angle %2d/%d] sigma_eff = %.5g s', j, N_angles, sigma_eff);
        if Lh < Lh_full
            T_full  = Lh_full/fs; T_kept = Lh/fs;
            energy_frac = erf(T_kept/(sqrt(2)*sigma_eff));
            fprintf(' | truncated %.3fs -> %.3fs (%.2f%% energy kept)', T_full, T_kept, 100*energy_frac);
            if energy_frac < 0.999, fprintf(' <-- WARNING: window far from Gaussian'); end
        end
        fprintf('\n');

        tau = (-Lh:Lh)' / fs;
        g_rot  = A_j * exp(-pi*b_j*tau.^2);
        dg_rot = (-2*pi*b_j) .* tau .* g_rot;

        offset = Lh_max - Lh;   % always >= 0 since Lh <= Lh_max
        g_eff(offset+1 : offset+2*Lh+1)  = g_eff(offset+1 : offset+2*Lh+1)  + g_rot;
        dg_eff(offset+1 : offset+2*Lh+1) = dg_eff(offset+1 : offset+2*Lh+1) + dg_rot;

        [STFT_j, ~, ~, ~] = compute_stft(x, g_rot, dg_rot, fs, nfft, hop, false);  % false = skip unused derivative FFT

        C_g_j = max(abs(sum(g_rot)), eps);
        mag_j = single(abs(STFT_j) / C_g_j);

        STFT_fused_sum = STFT_fused_sum + (mag_j.^2);
        n_included = n_included + 1;
        clear STFT_j mag_j;
    end

    fprintf('======================================================\n\n');
    STFT_fused = sqrt(STFT_fused_sum / n_included);

    info.N_angles = N_angles;
    info.shape_a  = a_eff;
    info.sigma_max = sigma_max;

    g_eff  = g_eff  / N_angles;   
    dg_eff = dg_eff / N_angles;

    [STFT_fused_coherent, dstft_fused_coherent, ~, ~] = ...
        compute_stft(x, g_eff, dg_eff, fs, nfft, hop, true);

    C_g_eff = abs(sum(g_eff));           
    g0_eff  = abs(g_eff(Lh_max + 1));    
    
end


%% MEMORY-EFFICIENT COMPUTE_STFT WITH HOPPING (derivative spectrogram now optional)
function [STFT, dstft, f, t] = compute_stft(x, g, dg, fs, nfft, hop, need_dstft)

    if nargin < 6 || isempty(hop), hop = 5; end
    if nargin < 7 || isempty(need_dstft), need_dstft = true; end

    x = x(:); g = g(:); dg = dg(:);
    N = length(x); Lw = length(g);
    if length(dg) ~= Lw, error('g and dg must have the same length.'); end
    if mod(Lw, 2) == 0, error('The window length must be odd.'); end

    Lh = floor(Lw / 2);
    if nfft < Lw, error('nfft must be greater than or equal to window length.'); end

    sample_indices = 1:hop:N;
    N_time = length(sample_indices);

    STFT  = zeros(nfft, N_time, 'single');
    dstft = zeros(nfft, N_time, 'single');
    x_pad = [zeros(Lh, 1); x; zeros(Lh, 1)];
    buf_g  = zeros(nfft, 1);
    buf_dg = zeros(nfft, 1);

    for idx = 1:N_time
        n = sample_indices(idx);
        segment = x_pad(n : n + 2*Lh);

        buf_g(:) = 0;
        buf_g(1:Lh+1) = segment(Lh+1:end) .* g(Lh+1:end);
        buf_g(end-Lh+1:end) = segment(1:Lh) .* g(1:Lh);
        STFT(:, idx) = single(fftshift(fft(buf_g)));

        if need_dstft
            buf_dg(:) = 0;
            buf_dg(1:Lh+1) = segment(Lh+1:end) .* dg(Lh+1:end);
            buf_dg(end-Lh+1:end) = segment(1:Lh) .* dg(1:Lh);
            dstft(:, idx) = single(fftshift(fft(buf_dg)));
        end
    end

    f = (-floor(nfft/2) : ceil(nfft/2)-1).' * fs / nfft;
    t = (sample_indices - 1) / fs;
end

%% HELPER FUNCTION FOR RENYI ENTROPY
function H = calculate_renyi_entropy(STFT_mag, alpha)
    P = STFT_mag.^2;
    total_energy = sum(P(:));
    if total_energy <= 0
        H = Inf;
        return;
    end
    P = P / total_energy;
    P = P(P > 0);

    if abs(alpha - 1) < 1e-12
        H = -sum(P .* log2(P));
    else
        H = log2(sum(P.^alpha)) / (1 - alpha);
    end
end

%% MODIFIED CONFIGURATION
function config = get_default_parameters()
    config.fs = 5000;
    config.file_name = 'faulty_stator_v00.mat';
    config.hop = 5;
    config.peak_thr = 0.01;
    config.peak_min_dist = 5;
    config.f_disp = 220;
    config.window_sigma_candidates = ...
        [0.005 0.01 0.015 0.02 0.025 0.03 0.04 0.05 0.055 0.06 0.065 0.07 0.075 0.080 0.085 0.09 0.095 0.1 0.105 0.11 0.115 0.120];
    config.window_opt_max_samples = 10000;
    config.window_opt_nfft = 4096;

    % Multi-Window Gauss Fusion Parameters
    config.multi_N_angles  = 12;      % number of orientations fused
    config.multi_shape_a   = 2.0;     % width RATIO between narrowest/widest orientation (1 = no diversity, degenerate; try 2-5)
    config.multi_sigma_max = 0.1;    % physical width (s) of the WIDEST-oriented window; keep within what nfft/fs can hold

    config.renyi_alpha = 2;
    config.plot_window_optimization = true;
    config.nfft = 4096;               % back to normal size now that width is controlled independently
    config.plot_signal = true;
    config.plot_fft = true;
    config.plot_stft = true;
    config.stft_amplitude_max = 1.2;
    
    config.signal_plot_xlim = [2 3];
end

%% plot_signal (Unchanged)
function plot_signal(t, x, config)
    figure;
    plot(t, real(x), 'LineWidth', 1.1);
    hold on;
    plot(t, imag(x), '--', 'LineWidth', 1.1);
    xlabel('Time (s)');
    ylabel('Amplitude');
    title('Complex Signal x(t)');
    legend('Real part', 'Imaginary part');
    grid on;
    if isfield(config, 'signal_plot_xlim') && ~isempty(config.signal_plot_xlim)
        xlim(config.signal_plot_xlim);
    end
end

%% detect_frequency_FFT (Unchanged)
function [f_pos, A_pos, f_neg, A_neg, fax, fft_orig] = detect_frequency_FFT(sig_x, config)
    N = length(sig_x);
    fft_orig = abs(fftshift(fft(sig_x)) / N);
    fax = (-floor(N/2):ceil(N/2)-1) * (config.fs/N);

    if isempty(config.fs)
        search_win = (fax ~= 0);
    else
        search_win = (fax >= -config.fs/2) & (fax <= config.fs/2) & (fax ~= 0);
    end

    f_search   = fax(search_win);
    mag_search = fft_orig(search_win);

    peak_threshold = config.peak_thr * max(fft_orig);
    [peaks, locs] = findpeaks(mag_search, f_search, ...
        'MinPeakHeight', peak_threshold, 'MinPeakDistance', config.peak_min_dist);

    f_pos = 0; A_pos = 0; f_neg = 0; A_neg = 0;

    pos_idx = locs > 0;
    if any(pos_idx)
        pos_peaks = peaks(pos_idx);
        pos_locs  = locs(pos_idx);
        [A_pos, idx] = max(pos_peaks);
        f_pos = pos_locs(idx);
    end

    neg_idx = locs < 0;
    if any(neg_idx)
        neg_peaks = peaks(neg_idx);
        neg_locs  = locs(neg_idx);
        [A_neg, idx] = max(neg_peaks);
        f_neg = neg_locs(idx);
    end
end

%% load_signal (Unchanged)
function [t, signal] = load_signal(config)
    data = load(config.file_name); 
    N = length(data.Te_f_12);
    fs = config.fs;
    t = (0:N-1) / fs;
    idx = (t >= 2);
    t = t(idx);

    ia = data.isa_f_12(idx, 1);
    ib = data.isb_f_12(idx, 1);
    ic = data.isc_f_12(idx, 1);

    signal = sqrt(2/3) * (ia + ib * exp(1i*2*pi/3) + ic * exp(1i*4*pi/3));
end

%% optimize_gaussian_window 
function [win, config_opt, result] = optimize_gaussian_window(sig_x, config)
    sig_x = sig_x(:);
    N = length(sig_x);
    fs = config.fs;

    if isfield(config, 'window_opt_max_samples') && N > config.window_opt_max_samples
        x_opt = sig_x(1:config.window_opt_max_samples);
    else
        x_opt = sig_x;
    end

    N_opt = length(x_opt);
    sigma_candidates = config.window_sigma_candidates;
    n_sigma = length(sigma_candidates);

    renyi_entropy = zeros(1, n_sigma);
    window_lengths = zeros(1, n_sigma);

    alpha = 2;
    if isfield(config, 'renyi_alpha'), alpha = config.renyi_alpha; end
    nfft = 2048;
    if isfield(config, 'window_opt_nfft'), nfft = config.window_opt_nfft; end

    for k = 1:n_sigma
        sigma = sigma_candidates(k);
        Lh = ceil(4 * sigma * fs);
        Lh_max_signal = floor((N_opt - 1) / 2);
        Lh = min(Lh, Lh_max_signal);

        if Lh < 1
            renyi_entropy(k) = Inf;
            continue;
        end

        winLen = 2 * Lh + 1;
        if winLen > nfft
            renyi_entropy(k) = Inf;
            continue;
        end

        tau = (-Lh:Lh)' / fs;
        g = exp(-0.5 * (tau / sigma).^2);
        dg = -tau / (sigma^2) .* g;

        window_lengths(k) = winLen;
        [STFT, ~, ~, ~] = compute_stft(x_opt, g, dg, fs, nfft, config.hop);

        P = abs(STFT).^2;
        total_energy = sum(P(:));
        if total_energy <= 0
            renyi_entropy(k) = Inf;
            continue;
        end

        P = P / total_energy;
        P = P(P > 0);

        if abs(alpha - 1) < 1e-12
            H = -sum(P .* log2(P));
        else
            H = log2(sum(P.^alpha)) / (1 - alpha);
        end
        renyi_entropy(k) = H;
    end

    [min_entropy, best_idx] = min(renyi_entropy);
    sigma_opt = sigma_candidates(best_idx);

    Lh_opt = ceil(4 * sigma_opt * fs);
    tau_opt = (-Lh_opt:Lh_opt)' / fs;
    g_opt = exp(-0.5 * (tau_opt / sigma_opt).^2);
    dg_opt = -tau_opt / (sigma_opt^2) .* g_opt;

    config_opt = config;
    config_opt.sigma = sigma_opt;
    config_opt.Lh = Lh_opt;
    config_opt.winLen = 2 * Lh_opt + 1;
    config_opt.gaussian_window = g_opt;
    config_opt.gaussian_derivative = dg_opt;

    result.sigma_candidates = sigma_candidates;
    result.renyi_entropy = renyi_entropy;
    result.window_lengths = window_lengths;
    result.best_idx = best_idx;
    result.sigma_opt = sigma_opt;
    result.min_entropy = min_entropy;

    fprintf('\nGaussian window optimization\n');
    fprintf('Number of candidates  = %d\n', n_sigma);
    fprintf('Optimal sigma          = %.6f s\n', sigma_opt);
    fprintf('Optimal window length = %d samples\n', config_opt.winLen);
    fprintf('Minimum Renyi entropy = %.6f\n', min_entropy);

    if isfield(config, 'plot_window_optimization') && config.plot_window_optimization
        figure;
        plot(sigma_candidates, renyi_entropy, 'o-', 'LineWidth', 1.2);
        hold on;
        plot(sigma_opt, min_entropy, 'rx', 'MarkerSize', 10, 'LineWidth', 2);
        xlabel('\sigma (seconds)');
        ylabel('Rényi entropy');
        title('Gaussian Window Optimization');
        grid on;
        legend('Candidate windows', 'Optimal window');
    end

    win = g_opt;
end

%% plot_fft (Unchanged)
function plot_fft(fax, fft_orig, config, f_pos, f_neg)
    figure;
    plot(fax, fft_orig, 'LineWidth', 1.2);
    hold on;
    if f_pos ~= 0
        A_pos = interp1(fax, fft_orig, f_pos, 'nearest');
        plot(f_pos, A_pos, 'o', 'MarkerSize', 8, 'LineWidth', 1.5);
    end
    if f_neg ~= 0
        A_neg = interp1(fax, fft_orig, f_neg, 'nearest');
        plot(f_neg, A_neg, 'o', 'MarkerSize', 8, 'LineWidth', 1.5);
    end
    xlabel('Frequency (Hz)');
    ylabel('Amplitude');
    title('Two-Sided FFT');
    grid on;
    xlim([-config.f_disp config.f_disp]);
end

%% plot_stft (Unchanged)
function plot_stft(t, f, STFT_phys, config, f_pos, f_neg)
    figure;
    imagesc(t, f, STFT_phys);
    axis xy;
    xlabel('Time (s)');
    ylabel('Frequency (Hz)');
    title('STFT Physical Amplitude');
    colorbar;
    ylim([-config.f_disp config.f_disp]);

    if isfield(config, 'stft_amplitude_max') && ~isempty(config.stft_amplitude_max)
        caxis([0 config.stft_amplitude_max]);
    end
    hold on;

    if f_pos ~= 0
        yline(f_pos, '--', 'LineWidth', 1.2);
    end
    if f_neg ~= 0
        yline(f_neg, '--', 'LineWidth', 1.2);
    end
end

%% SYNCHROSQUEEZED STFT (SST): FREQUENCY-AXIS REASSIGNMENT
% Input: a complex STFT and its matching derivative-window STFT (dstft).
% Reassigns each coefficient along frequency only:
%   omega_hat(t,eta) = eta - (1/2*pi) * Im( dstft(t,eta) / STFT(t,eta) )
function SSTFT_out = compute_sst(STFT, dstft, freqs, gamma_rel)

    if nargin < 4 || isempty(gamma_rel), gamma_rel = 1e-4; end

    [nfft, N_time] = size(STFT);
    freqs = freqs(:);
    df = freqs(2) - freqs(1);

    gamma = gamma_rel * max(abs(STFT(:)));
    valid = abs(STFT) >= gamma;

    freqMat   = repmat(freqs, 1, N_time);
    inst_freq = -Inf(nfft, N_time);
    ratio     = zeros(nfft, N_time);
    ratio(valid)     = dstft(valid) ./ STFT(valid);
    inst_freq(valid) = freqMat(valid) - imag(ratio(valid)) / (2*pi);

    in_range = valid & inst_freq >= freqs(1) & inst_freq <= freqs(end);
    target_k = zeros(nfft, N_time);
    target_k(in_range) = round((inst_freq(in_range) - freqs(1)) / df) + 1;
    target_k(in_range) = min(max(target_k(in_range), 1), nfft);

    SSTFT_out = zeros(nfft, N_time, 'like', STFT);
    for n = 1:N_time
        idx = in_range(:, n);
        if any(idx)
            SSTFT_out(:, n) = accumarray(target_k(idx, n), STFT(idx, n), [nfft, 1]);
        end
    end
end

%% plot_reconstruction (NEW)
function plot_reconstruction(t, x_ref, x_rec, ttl, ylims)
    figure;
    plot(t, real(x_ref), 'k', 'LineWidth', 1.2); hold on;
    plot(t, real(x_rec), 'r--', 'LineWidth', 1);
    xlabel('Time (s)'); ylabel('Real Amplitude');
    legend('Original (sampled)', 'Reconstructed');
    title(ttl);
    if nargin >= 5 && ~isempty(ylims)
        ylim(ylims);
    end
    xlim([2 3]);
    grid on;
end

%% RECONSTRUCTION METRICS (NEW)
function m = compute_recon_metrics(x_ref, x_rec)
    x_ref = x_ref(:); x_rec = x_rec(:);
    err = x_ref - x_rec;

    m.RMSE  = sqrt(mean(abs(err).^2));
    m.NRMSE = m.RMSE / sqrt(mean(abs(x_ref).^2));         % global L2-ratio (unitless)

    rel_mask = abs(x_ref) > 1e-6 * max(abs(x_ref));       % guard near-zero crossings
    m.RelErr_pct = mean(abs(err(rel_mask)) ./ abs(x_ref(rel_mask))) * 100;  % pointwise mean, %

    xr = x_ref - mean(x_ref);
    xc = x_rec - mean(x_rec);
    denom = norm(xr) * norm(xc);
    if denom > 0
        m.Corr = abs(sum(xr .* conj(xc))) / denom;        % complex correlation coefficient, [0,1]
    else
        m.Corr = NaN;
    end
end


function [t, x] = synthetic_signal()
    %% 1. Construction d'un signal complexe (chirp I/Q)
    fs = 5000;              % frequence d'echantillonnage (Hz)
    t  = (0:1/fs:2-1/fs)';  % duree 5 s

    % Chirp complexe : frequence instantanee variant de 50 Hz a 200 Hz
    f0 = 50; f1 = 200; T = t(end);
    k  = (f1 - f0) / T;                     % taux de balayage (Hz/s)
    phase = 2*pi*(f0*t + k/2*t.^2);
    x1 = exp(1i*phase);                     % composante 1 : chirp montant

    % Deuxieme composante : frequence NEGATIVE constante (impossible a
    % distinguer avec un signal reel, mais naturel avec un signal complexe)
    x2 = 0.6*exp(-1i*2*pi*80*t);

    x = x1 + x2;                            % signal complexe final (bruit ajoute)
    x = x + 0.02*(randn(size(x)) + 1i*randn(size(x)));
end