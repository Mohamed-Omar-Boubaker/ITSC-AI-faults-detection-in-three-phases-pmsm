%% MAIN FUNCTION
function Custom_SWT_energy_extraction()
    t_total = tic;
    
    % 1. Configuration & Data Loading
    p = params();
    
    % Choose testing signal (uncomment the desired one)
    [t, x] = load_faulty_signal(p);         % Faulty sample from DB ~±33Hz
    %[t, x] = load_healthy_signal(p);        % Healthy sample from DB ~+33Hz
    %[t, x, IF] = load_signal_freq();        % Signal with transient freq and IF
    
    if p.resample
        % Resample signal to lower frequency        
        [t, x, p.fs] = resample_signal(t, x, p.fs, p.fs_resample);
    end
    
    IF = -200 * ones(length(t), 1);
    p.N = length(x);
    
    % Wavelet scale and grid setup
    [scales, f] = morseScales(p.N, p.fs, p.beta_opt, p.gamma, p.voices_per_octave);
    [da, ~] = scaleSpacing(scales, f);
    
    % 2. Bilateral Continuous Wavelet Transform (CWT)
    fprintf("Running Native 3D Bilateral CWT...\n");
    t_step = tic;
    [W, dW] = compute_CWT(x, p, scales);
    fprintf("  -> CWT took %.3f seconds\n", toc(t_step));
    
    if p.plotting
        plot_CWT_scalogram(W, f, p, t);
    end
    
    % 3. Synchrosqueezed Wavelet Transform (SWT)
    fprintf("Running Synchrosqueezing...\n");
    t_step = tic;
    [Tx_pos, f] = SWT(W(:,:,1), dW(:,:,1), f, scales, p.thresh_ratio, da);
    [Tx_neg, f] = SWT(W(:,:,2), dW(:,:,2), f, scales, p.thresh_ratio, da);
    fprintf("  -> SWT took %.3f seconds\n", toc(t_step));
    
    clear dW; % Free CWT derivative matrix
    
    % Morse wavelet admissibility constant
    C_psi = admissibilityConstant(p.beta_opt, p.gamma);
    
    if p.plot_SWT
        plot_SWT_scalogram(Tx_pos, Tx_neg, f, C_psi, t, p);
    end
    
    % 4. Signal Reconstruction Check (Optional)
    if p.reconstruct_original_signal
        x_hat_cwt = (sum(W(:,:,1) .* (da./scales), 1) + sum(W(:,:,2) .* (da./scales), 1)) / C_psi;
        x_hat_swt = (sum(Tx_pos, 1) + sum(Tx_neg, 1)) / C_psi;
        
        plot_reconstructed_signals(t, x, x_hat_cwt, x_hat_swt);
        
        if p.show_metrics
            compute_signal_metrics(x, x_hat_cwt);
            compute_signal_metrics(x, x_hat_swt);
        end
        
        clear x_hat_cwt x_hat_swt;
    end
    
    clear W; % Free CWT coefficient matrix
    
    % 5. Build Combined Bilateral Matrices
    [Tx_all, f_axis, Tx_all_complex] = build_bilateral_matrices(Tx_pos, Tx_neg, f, C_psi);
    clear Tx_pos Tx_neg;
    
    % 6. Dual-Limit Squeeze Energy Extraction
    fprintf("Running Dual-Limit Squeeze...\n");
    t_step = tic;
    [rec_energy, target_times, WE_cyclic, ratio_cyclic, WE_global, ratio_global] = dual_limit_squeeze_energy(...
        Tx_all, Tx_all_complex, f_axis, t, p, IF, C_psi);
    fprintf("  -> Energy extraction took %.3f seconds\n", toc(t_step));
    
    % 7. Display Extracted Features
    fprintf("\n========== Extracted Features ==========\n");
    fprintf("Global Wavelet Entropy        : %.6f bits\n", WE_global);
    fprintf("Tracked Band Energy Ratio     : %.6f\n", ratio_global);
    fprintf("Mean Cyclic Energy            : %.6f\n", mean(rec_energy));
    fprintf("Mean Cyclic Wavelet Entropy   : %.6f bits\n", mean(WE_cyclic));
    fprintf("Mean Cyclic Energy Ratio      : %.6f\n", mean(ratio_cyclic));
    fprintf("=======================================\n");
    
    fprintf("\nDone! Total execution time: %.3f seconds\n", toc(t_total));
end

%% CONFIGURATION PARAMETERS
function p = params()
    p.fs                             = 100000;
    p.fs_resample                    = 10000;
    p.N                              = 0;
    p.gamma                          = 3;
    p.beta_opt                       = 190;
    p.voices_per_octave              = 14;
    p.f_disp                         = 500; 
    p.thresh_ratio                   = 1e-4;
    p.resample                       = true;
    p.plotting                       = false;
    p.plot_SWT                       = true;
    p.plot_rec_sig                   = false;
    p.show_metrics                   = true;
    p.reconstruct_original_signal    = false;
    p.plot_hard_threshold_SWT_direct = false;
    p.plot_hard_threshold_SWT_rec    = false;
    p.plot_squeeze_SWT_direct        = false;   
    p.plot_squeeze_SWT_rec           = false;   
    p.plot_cyclic_energy             = true;
    p.plot_cyclic_wavelet_entropy    = true;
    p.plot_cyclic_energy_ratio       = true;
    p.search_radius_hz               = 5;
    p.max_expansion_hz               = 3;
    p.window_sec                     = 0.2;
end

%% CWT COMPUTATION
function [W, dW] = compute_CWT(x, p, scales)
    x = x(:).';
    Fs = p.fs;
    N = length(x);
    
    X = fft(x);
    w = (2*pi/N) * (0:N-1);
    omega_rad = (2*pi*Fs/N) * [0:floor((N-1)/2), -floor(N/2):-1];
    
    Ns = length(scales);
    W  = zeros(Ns, N, 2);
    dW = zeros(Ns, N, 2);
    
    for k = 1:Ns
        a = scales(k);
        % Positive frequencies
        psi_pos = morseWaveletFFT(w, a, p.beta_opt, p.gamma, 'pos');
        W(k,:,1)  = ifft(X .* conj(psi_pos));
        dW(k,:,1) = ifft(X .* (1i * omega_rad) .* conj(psi_pos));
        
        % Negative frequencies
        psi_neg = morseWaveletFFT(w, a, p.beta_opt, p.gamma, 'neg');
        W(k,:,2)  = ifft(X .* conj(psi_neg));
        dW(k,:,2) = ifft(X .* (1i * omega_rad) .* conj(psi_neg));
    end
end

%% MORSE WAVELET HELPER FUNCTIONS
function [scales, f] = morseScales(N, Fs, beta, gamma, VPO)
    w0 = (beta / gamma)^(1 / gamma);
    fmax = Fs / 2;
    fmin = Fs / N;
    octaves = log2(fmax / fmin);
    Nv = floor(octaves * VPO);
    
    f = fmax ./ 2.^((0:Nv)/VPO);
    f = f(f >= fmin).';
    scales = (w0 * Fs) ./ (2 * pi * f);
end

function psi = morseWaveletFFT(w, scale, beta, gamma, side)
    w_s = scale .* w;
    if strcmp(side, 'pos')
        w_domain = max(w_s, 0);
        w_domain(w > pi) = 0;
    else
        w_domain = zeros(size(w));
        neg_idx = w > pi;
        w_domain(neg_idx) = scale * (2*pi - w(neg_idx));
    end
    
    psi = zeros(size(w_domain));
    idx = w_domain > 0;
    
    logA = (beta/gamma) * (1 + log(gamma/beta));
    log_psi0 = logA + beta .* log(w_domain(idx)) - (w_domain(idx)).^gamma;
    psi(idx) = exp(log_psi0);
end

function C_psi = admissibilityConstant(beta, gamma) 
    logA = (beta/gamma) * (1 + log(gamma/beta));      
    logC = logA - log(gamma) + gammaln(beta/gamma);
    C_psi = exp(logC);
end

function [da, dw] = scaleSpacing(scales, f)
    scales = scales(:); 
    f = f(:);
    w_rad = 2*pi*f;
    Ns = length(scales);

    da = zeros(Ns, 1);
    da(1)       = abs(scales(2) - scales(1));
    da(2:end-1) = abs(scales(3:end) - scales(1:end-2)) / 2;
    da(end)     = abs(scales(end) - scales(end-1));

    dw = zeros(Ns, 1);
    dw(1)       = abs(w_rad(2) - w_rad(1));
    dw(2:end-1) = abs(w_rad(3:end) - w_rad(1:end-2)) / 2;
    dw(end)     = abs(w_rad(end) - w_rad(end-1));
end

%% SYNCHROSQUEEZED WAVELET TRANSFORM (SWT)
function [Tx, f_axis] = SWT(W, dW, f, scales, thresh_ratio, da)
    if nargin < 5 || isempty(thresh_ratio), thresh_ratio = 0.05; end
    f = f(:); 
    scales = scales(:);
    [Ns, N] = size(W);

    w_rad = 2 * pi * f;
    [w_sorted, sort_idx] = sort(w_rad, 'ascend');
    row_lookup = (1:Ns)';

    thresh = thresh_ratio * max(abs(W), [], 'all');
    valid = abs(W) > thresh;
    valid_idx = find(valid);

    if isempty(valid_idx)
        Tx = zeros(Ns, N, 'like', W);
        f_axis = f;
        return;
    end

    % Instantaneous frequency estimation
    IF_rad = zeros(Ns, N);
    IF_rad(valid) = abs(imag(dW(valid) ./ W(valid)));

    % Interpolation to scale bins
    rows = interp1(w_sorted, row_lookup(sort_idx), IF_rad(valid_idx), 'nearest', 'extrap');
    rows = max(1, min(Ns, rows));

    % Integration weighting
    scales_mat = repmat(scales, 1, N);
    da_mat     = repmat(da, 1, N);
    weights    = W(valid_idx) .* (scales_mat(valid_idx) .^ (-1)) .* da_mat(valid_idx);

    [~, cols_idx] = ind2sub([Ns, N], valid_idx);
    Tx = accumarray([rows, cols_idx], weights, [Ns, N]);
    f_axis = f;
end

function [Tx_all, f_axis, Tx_all_complex] = build_bilateral_matrices(Tx_pos, Tx_neg, f, C_psi)
    NumTime = size(Tx_pos, 2);
    f_axis = [-f; 0; flipud(f)];
    [f_axis, sortIdx] = sort(f_axis, 'ascend');

    % Magnitudes (visualization)
    Tx_all = [abs(Tx_neg); zeros(1, NumTime); flipud(abs(Tx_pos))] / C_psi;
    Tx_all = Tx_all(sortIdx, :);

    % Complex (reconstruction)
    Tx_all_complex = [Tx_neg; zeros(1, NumTime); flipud(Tx_pos)];
    Tx_all_complex = Tx_all_complex(sortIdx, :);
end

%% ADAPTIVE SQUEEZE & ENERGY EXTRACTION
function [rec_energy, cyclic_time, WE_cyclic, ratio_cyclic, WE_global, ratio_global] = dual_limit_squeeze_energy(...
    Tx_all, Tx_all_complex, f_axis, t, p, IF, C_psi)

    % 1. Adaptive tracking on SWT domain
    [band_min_f, band_max_f_vis, instant_energy, band_indices, lower_thresh_history] = ...
        track_adaptive_squeeze_bands(Tx_all, f_axis, IF, p.search_radius_hz, p.max_expansion_hz);

    [WE_global, WE_cyclic, cyclic_time] = compute_wavelet_entropy(Tx_all, f_axis, t, p);
    [ratio_cyclic, ratio_global] = compute_subband_energy_ratio(Tx_all, band_indices, lower_thresh_history, t, p);

    % Energy extraction via signal reconstruction
    [x_rec, instant_energy_rec] = reconstruct_squeeze_signal_from_band(...
        Tx_all_complex, Tx_all, band_indices, lower_thresh_history, C_psi);

    clear Tx_all_complex; % Free complex spectrum matrix

    [cyclic_time, rec_energy] = compute_cyclic_energy(instant_energy_rec, t, p.window_sec);

    if p.plot_squeeze_SWT_rec
        plot_squeeze_reconstruction_results(t, f_axis, IF, Tx_all, band_min_f, band_max_f_vis, ...
            x_rec, instant_energy_rec, cyclic_time, rec_energy, p);
    end
    
    if isfield(p, 'plot_cyclic_energy') && p.plot_cyclic_energy
        plot_cyclic_energy_fig(cyclic_time, rec_energy);
    end
end

function [band_min_f, band_max_f_vis, instant_energy, band_indices, lower_thresh_history] = ...
    track_adaptive_squeeze_bands(Tx_all, f_axis, IF, search_radius_hz, max_expansion_hz)

    [NumFreq, NumTime] = size(Tx_all);

    band_min_f           = zeros(NumTime, 1);
    band_max_f_vis       = zeros(NumTime, 1);
    instant_energy       = zeros(NumTime, 1);
    band_indices         = cell(NumTime, 1);
    lower_thresh_history = zeros(NumTime, 1);

    shrinkage_factor = 0.125;           
    max_iters        = 100;     

    for idx = 1:NumTime
        current_IF = IF(idx);
        [~, nominal_f_idx] = min(abs(f_axis - current_IF));

        search_range = find(f_axis >= (current_IF - search_radius_hz) & ...
                            f_axis <= (current_IF + search_radius_hz));

        if isempty(search_range)
            search_range = nominal_f_idx;
        end
        slice_amps = Tx_all(search_range, idx);

        [peak_amp, local_peak_offset] = max(slice_amps);
        center_f_idx = search_range(local_peak_offset);
        f_center     = f_axis(center_f_idx);

        if peak_amp < eps
            band_min_f(idx)     = f_center;
            band_max_f_vis(idx) = f_axis(min(NumFreq, center_f_idx + 1));
            band_indices{idx}   = center_f_idx;
            continue;
        end

        upper_thresh = 0.50 * peak_amp;
        lower_thresh = 0.01 * peak_amp;

        for iter = 1:max_iters
            N_uth = count_peaks_above_thresh(slice_amps, upper_thresh);
            N_lth = count_peaks_above_thresh(slice_amps, lower_thresh);

            if N_uth == N_lth || (upper_thresh - lower_thresh) < 1e-7 * peak_amp
                break;
            end
            
            D = upper_thresh - lower_thresh;
            upper_thresh = upper_thresh - (shrinkage_factor * D);
            lower_thresh = lower_thresh + (shrinkage_factor * D);
        end

        lower_thresh_history(idx) = lower_thresh;
        
        limit_up   = find(f_axis <= (current_IF + max_expansion_hz), 1, 'last');
        limit_down = find(f_axis >= (current_IF - max_expansion_hz), 1, 'first');

        if isempty(limit_up),   limit_up   = NumFreq; end
        if isempty(limit_down), limit_down = 1;       end

        up_idx = center_f_idx;
        while up_idx < limit_up && Tx_all(up_idx + 1, idx) >= lower_thresh
            up_idx = up_idx + 1;
        end

        down_idx = center_f_idx;
        while down_idx > limit_down && Tx_all(down_idx - 1, idx) >= lower_thresh
            down_idx = down_idx - 1;
        end

        down_idx = max(1, down_idx);
        up_idx   = min(NumFreq, up_idx);

        band_min_f(idx) = f_axis(down_idx);
        if up_idx < NumFreq
            band_max_f_vis(idx) = f_axis(up_idx + 1);
        else
            band_max_f_vis(idx) = f_axis(up_idx);
        end

        band_range        = down_idx:up_idx;
        band_indices{idx} = band_range;

        valid_coefs = Tx_all(band_range, idx);
        valid_coefs(valid_coefs < lower_thresh) = 0;
        instant_energy(idx) = sum(valid_coefs.^2);
    end
end

function num_peaks = count_peaks_above_thresh(signal_slice, thresh)
    above_idx = find(signal_slice > thresh);
    if isempty(above_idx)
        num_peaks = 0;
        return;
    end

    len = length(signal_slice);
    if len == 1
        num_peaks = 1;
        return;
    end

    is_peak = false(size(signal_slice));
    for i = 2:(len - 1)
        if signal_slice(i) > thresh && ...
           signal_slice(i) >= signal_slice(i-1) && ...
           signal_slice(i) >= signal_slice(i+1)
            is_peak(i) = true;
        end
    end

    if signal_slice(1) > thresh && signal_slice(1) > signal_slice(2)
        is_peak(1) = true;
    end
    if signal_slice(end) > thresh && signal_slice(end) > signal_slice(end-1)
        is_peak(end) = true;
    end

    num_peaks = sum(is_peak);
    if num_peaks == 0 && ~isempty(above_idx)
        num_peaks = 1;
    end
end

function [x_rec, instant_energy_rec] = reconstruct_squeeze_signal_from_band(...
    Tx_all_complex, Tx_all, band_indices, lower_thresh_history, C_psi)

    NumTime = length(band_indices);
    x_rec = zeros(NumTime, 1);
    instant_energy_rec = zeros(NumTime, 1);

    for idx = 1:NumTime
        band_range = band_indices{idx};
        lower_thresh = lower_thresh_history(idx);

        coefs_tx_complex = Tx_all_complex(band_range, idx);
        valid_swt_check  = Tx_all(band_range, idx);

        coefs_tx_complex(valid_swt_check < lower_thresh) = 0;

        x_rec(idx) = sum(coefs_tx_complex, 1) / C_psi;
        instant_energy_rec(idx) = abs(x_rec(idx))^2;
    end
end

function [cyclic_time, cyclic_energy] = compute_cyclic_energy(instant_energy, t, window_sec)
    dt = t(2) - t(1);
    N = length(t);
    samples_per_window = round(window_sec / dt);
    num_cycles = ceil(N / samples_per_window);
    
    cyclic_energy = zeros(num_cycles, 1);
    cyclic_time   = zeros(num_cycles, 1);

    for c = 1:num_cycles
        start_samp = (c-1) * samples_per_window + 1;
        end_samp   = min(c * samples_per_window, N);
        cyclic_energy(c) = sum(instant_energy(start_samp:end_samp));
        cyclic_time(c)   = mean(t(start_samp:end_samp));
    end
end

function [cyclic_vals, target_times] = track_cyclic_energy(instant_energy, t, p, label_tag)
    target_times = (t(1) : p.window_sec : t(end)).';
    num_points = length(target_times);
    cyclic_vals = zeros(num_points, 1);
    
    for i = 1:num_points
        [~, exact_idx] = min(abs(t - target_times(i)));
        cyclic_vals(i) = instant_energy(exact_idx);
    end
    
    if isfield(p, 'plot_cyclic_energy') && p.plot_cyclic_energy
        plot_cyclic_energy_fig(target_times, cyclic_vals);
    end
end

%% DATA LOADING HELPERS
function [t, x] = load_faulty_signal(p)
    ia = double(h5read('pmsm_itsc_dataset.h5', '/current_phase_A', [1 100], [200000 1]));
    ib = double(h5read('pmsm_itsc_dataset.h5', '/current_phase_B', [1 100], [200000 1]));
    ic = double(h5read('pmsm_itsc_dataset.h5', '/current_phase_C', [1 100], [200000 1]));
    
    x = sqrt(2/3) * (ia + ib * exp(1i*2*pi/3) + ic * exp(1i*4*pi/3));
    
    N = length(ia);
    t = (0:N-1) / p.fs;

end

function [t, x] = load_healthy_signal(p)
    data = load("healthy_v00.mat"); 
    N = length(data.Te_h);
    t = (0:N-1) / p.fs;
    idx = (t >= 2);
    t = t(idx);

    ia = double(data.isa_h(idx, 1));
    ib = double(data.isb_h(idx, 1));
    ic = double(data.isc_h(idx, 1));

    x = sqrt(2/3) * (ia + ib * exp(1i*2*pi/3) + ic * exp(1i*4*pi/3));
end

function [t, x, IF] = load_signal_freq()
    fs      = 5000;     
    T_total = 4;          
    f0      = 20;         
    f1      = 30;         
    T_const = 1;          
    T_trans = T_total - 2*T_const;  
    A       = 0.3;        

    k = 10 / T_trans;
    t = (0:1/fs:T_total-1/fs).';   
    IF = zeros(size(t));

    idx1 = t < T_const;
    IF(idx1) = f0;

    idx2 = (t >= T_const) & (t < T_const + T_trans);
    t_rel = t(idx2) - T_const - T_trans/2;   
    sigmoid = 1 ./ (1 + exp(-k * t_rel));    
    IF(idx2) = f0 + (f1 - f0) * sigmoid;

    idx3 = t >= T_const + T_trans;
    IF(idx3) = f1;

    phase = 2*pi*cumtrapz(t, IF);
    x = A * sin(phase);
    x = hilbert(x);
end

%% METRICS & VISUALIZATIONS
function metrics = compute_signal_metrics(x, x_hat)
    x = x(:).';
    x_hat = x_hat(:).';
    err = x - x_hat;
    
    metrics.RMSE = sqrt(mean(abs(err).^2));
    metrics.NRMSE = norm(err) / norm(x);

    signal_energy = sum(abs(x).^2);
    error_energy  = sum(abs(err).^2);

    metrics.SNR = 10*log10(signal_energy/error_energy);
    metrics.Energy_original = signal_energy;
    metrics.Energy_reconstructed = sum(abs(x_hat).^2);

    p_dist = abs(x).^2;
    p_dist = p_dist / sum(p_dist);
    p_dist(p_dist == 0) = [];

    metrics.Shannon_entropy = -sum(p_dist .* log2(p_dist));

    alpha = 2;
    metrics.Renyi_entropy = 1/(1-alpha) * log2(sum(p_dist.^alpha));

    C = corrcoef(real(x), real(x_hat));
    metrics.Correlation = C(1,2);

    fprintf('\n========== Signal Metrics ==========\n');
    fprintf('RMSE                  : %.6e\n', metrics.RMSE);
    fprintf('NRMSE                 : %.6e\n', metrics.NRMSE);
    fprintf('SNR                   : %.3f dB\n', metrics.SNR);
    fprintf('Energy original       : %.6e\n', metrics.Energy_original);
    fprintf('Energy reconstructed  : %.6e\n', metrics.Energy_reconstructed);
    fprintf('Shannon entropy       : %.6f bits\n', metrics.Shannon_entropy);
    fprintf('Renyi entropy         : %.6f bits\n', metrics.Renyi_entropy);
    fprintf('Correlation           : %.12f\n', metrics.Correlation);
    fprintf('====================================\n');
end

function plot_CWT_scalogram(W, f, p, t)
    f_axis = [-f; 0; flipud(f)];
    amp_all = [abs(W(:,:,2)); zeros(1, p.N); flipud(abs(W(:,:,1)))];
    
    [f_axis, sortIdx] = sort(f_axis, 'ascend');
    amp_all = amp_all(sortIdx, :);

    keep = (f_axis >= -p.f_disp) & (f_axis <= p.f_disp);

    figure;
    surf(t, f_axis(keep), amp_all(keep, :), 'EdgeColor', 'none');
    view(2); axis xy;
    xlabel("Time (s)"); ylabel("Frequency (Hz)");
    title("Bilateral Wavelet Scalogram");
    ylim([-p.f_disp, p.f_disp]); 
    caxis([0 max(amp_all(:))]); colorbar;
    yline(0, '--w', 'LineWidth', 1.5);
end

function plot_SWT_scalogram(Tx_pos, Tx_neg, f, C_psi, t, p) 
    f_axis = [-f; 0; flipud(f)];
    [f_axis, sortIdx] = sort(f_axis, 'ascend');
    keep = (f_axis >= -p.f_disp) & (f_axis <= p.f_disp); 
    
    amp_all_swt = [abs(Tx_neg); zeros(1, p.N); flipud(abs(Tx_pos))] / C_psi;
    amp_all_swt = amp_all_swt(sortIdx, :);

    figure;
    surf(t, f_axis(keep), amp_all_swt(keep, :), 'EdgeColor', 'none');
    view(2); axis xy;
    xlabel("Time (s)"); ylabel("Frequency (Hz)");
    title("Bilateral Synchrosqueezed Wavelet Transform (SWT)");
    ylim([-p.f_disp, p.f_disp]);
    caxis([0 max(amp_all_swt(:))]); 
    colorbar;
    yline(0, '--w', 'LineWidth', 1.5);
end

function plot_reconstructed_signals(t, x, x_hat_cwt, x_hat_swt)   
    figure; 
    
    ax(1) = subplot(3, 1, 1);
    plot(t, real(x), 'k', 'LineWidth', 1.4);
    ylabel('Re\{x(t)\}');
    title('Original Signal');
    grid on;
    
    ax(2) = subplot(3, 1, 2);
    plot(t, real(x_hat_cwt), 'r', 'LineWidth', 1.2);
    ylabel('Re\{x(t)\}');
    title('Reconstructed Signal from CWT');
    grid on;
    
    ax(3) = subplot(3, 1, 3);
    plot(t, real(x_hat_swt), 'b', 'LineWidth', 1.2);
    xlabel('Time (s)');
    ylabel('Re\{x(t)\}');
    title('Reconstructed Signal from SWT');
    grid on;

    linkaxes(ax, 'xy');
    xlim([2 2.5]);
end

function plot_squeeze_reconstruction_results(t, f_axis, IF, Tx_all, band_min_f, band_max_f_vis, ...
    x_rec, instant_energy_rec, cyclic_time, cyclic_energy_rec, p)

    keep = (f_axis >= -p.f_disp) & (f_axis <= p.f_disp);

    figure('Color', 'w', 'Name', 'SWT Scalogram & Squeeze Band');
    surf(t, f_axis(keep), Tx_all(keep, :), 'EdgeColor', 'none');
    view(2); axis xy; hold on;

    z_top = ones(size(t)) * max(Tx_all(:)) * 2;
    h_if   = plot3(t, IF, z_top, 'w--', 'LineWidth', 2);
    h_band = plot3(t, band_min_f, z_top, 'r-', 'LineWidth', 1.5);
             plot3(t, band_max_f_vis, z_top, 'r-', 'LineWidth', 1.5);

    xlabel('Time (s)'); ylabel('Frequency (Hz)');
    title('SWT Scalogram with Selected Adaptive Band & SWT Inversion');
    ylim([-p.f_disp, p.f_disp]); colorbar; caxis([0 max(Tx_all(:))]);
    legend([h_if, h_band], {'True IF Profile', 'Adaptive Band Boundaries'}, 'Location', 'Best');
    hold off;

    figure('Color', 'w', 'Name', 'Reconstructed Signal & Temporal Energy');
    subplot(3,1,1);
    plot(t, real(x_rec), 'g', 'LineWidth', 1.2); grid on;
    xlabel('Time (s)'); ylabel('Amplitude');
    title('Denoised Reconstructed Signal (Real Part)'); xlim([2 2.5]);

    subplot(3,1,2);
    plot(t, instant_energy_rec, 'b', 'LineWidth', 1.5); grid on;
    xlabel('Time (s)'); ylabel('Instantaneous Energy');
    title('Instantaneous Temporal Energy inside Tracked Band');

    subplot(3,1,3);
    stem(cyclic_time, cyclic_energy_rec, 'r', 'Filled', 'LineWidth', 1.5); grid on;
    xlabel('Time (s)'); ylabel('Accumulated Energy');
    title('Cyclic Energy Profile of Reconstructed Signal');
end

%% SIGNAL RESAMPLING
function [t_rs, x_rs, fs_out] = resample_signal(t, x, fs_orig, fs_target)

    % Convert input to double precision (required by resample)
    x = double(x);
    t = double(t);

    % Calculate rational sample rate conversion factors
    [p_factor, q_factor] = rat(fs_target / fs_orig);

    % Resample signal with anti-aliasing filtering
    x_rs = resample(x, p_factor, q_factor);

    % Calculate effective target frequency and build matching time vector
    fs_out = fs_orig * (p_factor / q_factor);
    N_rs = length(x_rs);
    t_rs = t(1) + (0:N_rs-1)' / fs_out;

    % Preserve input vector orientation
    if isrow(x), x_rs = x_rs.'; end
    if isrow(t), t_rs = t_rs.'; end
end

%% WAVELET ENTROPY (SHANNON) — SWT-BASED FEATURE
function [entropy_global, entropy_cyclic, target_times] = compute_wavelet_entropy(Tx_all, f_axis, t, p)
   
    E_freq = sum(Tx_all.^2, 2);           % Energy per frequency bin
    entropy_global = shannon_entropy_from_energy(E_freq);

    % Time-resolved cyclic entropy over windows
    dt = t(2) - t(1);
    N  = length(t);
    samples_per_window = round(p.window_sec / dt);
    num_windows = ceil(N / samples_per_window);

    entropy_cyclic = zeros(num_windows, 1);
    target_times   = zeros(num_windows, 1);

    for w = 1:num_windows
        start_samp = (w-1) * samples_per_window + 1;
        end_samp   = min(w * samples_per_window, N);

        E_win = sum(Tx_all(:, start_samp:end_samp).^2, 2);
        entropy_cyclic(w) = shannon_entropy_from_energy(E_win);
        target_times(w)   = mean(t(start_samp:end_samp));
    end

    if isfield(p, 'plot_cyclic_wavelet_entropy') && p.plot_cyclic_wavelet_entropy
        plot_cyclic_wavelet_entropy(target_times, entropy_cyclic);
    end
end

function H = shannon_entropy_from_energy(E)
    E = E(:);
    total = sum(E);
    if total <= 0
        H = 0;
        return;
    end
    p_dist = E / total;
    p_dist(p_dist == 0) = [];
    H = -sum(p_dist .* log2(p_dist));
end

%% SUB-BAND ENERGY RATIO — SWT-BASED FEATURE
function [ratio_cyclic, ratio_global] = compute_subband_energy_ratio(Tx_all, band_indices, lower_thresh_history, t, p)

    NumTime = size(Tx_all, 2);
    band_energy_time  = zeros(NumTime, 1);
    total_energy_time = zeros(NumTime, 1);

    apply_mask = nargin >= 3 && ~isempty(lower_thresh_history);

    for idx = 1:NumTime
        band_range = band_indices{idx};
        coefs = Tx_all(band_range, idx);

        if apply_mask
            coefs(coefs < lower_thresh_history(idx)) = 0;
        end

        band_energy_time(idx)  = sum(coefs.^2);
        total_energy_time(idx) = sum(Tx_all(:, idx).^2);
    end

    ratio_global = sum(band_energy_time) / max(sum(total_energy_time), eps);

    % Cyclic Sub-Band Energy Ratio calculation over windows
    if nargin >= 5 && ~isempty(t) && ~isempty(p)
        dt = t(2) - t(1);
        N  = length(t);
        samples_per_window = round(p.window_sec / dt);
        num_windows = ceil(N / samples_per_window);

        ratio_cyclic = zeros(num_windows, 1);
        target_times = zeros(num_windows, 1);

        for w = 1:num_windows
            start_samp = (w-1) * samples_per_window + 1;
            end_samp   = min(w * samples_per_window, N);

            band_E_win  = sum(band_energy_time(start_samp:end_samp));
            total_E_win = sum(total_energy_time(start_samp:end_samp));

            if total_E_win > 0
                ratio_cyclic(w) = band_E_win / total_E_win;
            else
                ratio_cyclic(w) = 0;
            end
            target_times(w) = mean(t(start_samp:end_samp));
        end

        if isfield(p, 'plot_cyclic_energy_ratio') && p.plot_cyclic_energy_ratio
            plot_cyclic_energy_ratio(target_times, ratio_cyclic);
        end
    else
        ratio_cyclic = band_energy_time ./ max(total_energy_time, eps);
    end
end

%% FEATURE PLOTTING FUNCTIONS
function plot_cyclic_wavelet_entropy(target_times, WE_cyclic)
    figure('Color', 'w', 'Name', 'Cyclic Wavelet Entropy');
    plot(target_times, WE_cyclic, 'b-o', 'LineWidth', 1.5, 'MarkerFaceColor', 'b');
    grid on;
    xlabel('Time (s)');
    ylabel('Entropy (bits)');
    title('Cyclic Wavelet Entropy');
end

function plot_cyclic_energy_ratio(target_times, ratio_cyclic)
    figure('Color', 'w', 'Name', 'Cyclic Tracked Band Energy Ratio');
    plot(target_times, ratio_cyclic, 'r-s', 'LineWidth', 1.5, 'MarkerFaceColor', 'r');
    grid on;
    xlabel('Time (s)');
    ylabel('Energy Ratio');
    title('Cyclic Tracked Band Energy Ratio');
end

function plot_cyclic_energy_fig(target_times, rec_energy)
    figure('Color', 'w', 'Name', 'Cyclic Energy');
    plot(target_times, rec_energy, 'k-d', 'LineWidth', 1.5, 'MarkerFaceColor', 'k');
    grid on;
    xlabel('Time (s)');
    ylabel('Cyclic Energy');
    title('Mean / Cyclic Energy of Reconstructed Signal');
end