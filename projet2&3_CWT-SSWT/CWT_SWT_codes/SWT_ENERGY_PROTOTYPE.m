%% main 
function SWT_ENERGY_PROTOTYPE()
    close all;
    p = params();       % loading params

    % Define all candidate signals to process sequentially
    signal_list = {
        'Synthetic Signal (Pure 50Hz + -80Hz)', ...
        'Chirp Signal (Chirp 50->200Hz + -80Hz)', ...
        'Faulty Signal (Faulty sample from DB ~±33Hz)', ...
        'Healthy Signal (Healthy sample from DB ~+33Hz)', ...
        'Signal with Transient Frequency'
    };

    num_signals = length(signal_list);

    for s = 1:num_signals
        sig_name = signal_list{s};

        % Mention signal name clearly before starting treatment
        fprintf('\n========================================================\n');
        fprintf('  STARTING PROCESSING FOR SIGNAL %d/%d: %s\n', s, num_signals, sig_name);
        fprintf('========================================================\n\n');

        % Reset optional return variables for each signal run
        currents = [];
        IF = [];

        % Load selected signal based on iteration step
        switch s
            case 1
                [t, x] = synthetic_signal(p);
            case 2
                [t, x] = load_chirp_signal(p);
            case 3
                [t, x, currents] = load_faulty_signal(p);
            case 4
                [t, x, currents] = load_healthy_signal(p);
            case 5
                [t, x, IF] = load_signal_freq();
        end

        % Plot the three phase currents if present
        if ~isempty(currents)
            plot_three_phase_currents(t, currents, sprintf('[%s] Three-Phase Currents (Raw)', sig_name));

            % Save three-phase currents to a signal-unique CSV file
            csv_filename = sprintf('three_phase_currents_sig%d.csv', s);
            writematrix([t(:), currents.ia(:), currents.ib(:), currents.ic(:)], csv_filename);
            plot_csv_data(csv_filename, sprintf('[%s] Three-Phase Currents (from CSV)', sig_name));
        end

        % Add noise with a target SNR of 15 dB (lower dB = more noise)
        % x = add_custom_wgn(x, 15);

        % Set fallback frequency of interest IF if not assigned by loader
        if isempty(IF)
            IF = -33 * ones(length(t), 1);
        end

        % Signal length
        p.N = length(x);

        % FFT analysis and spectrum plot
        plot_FFT_spectrum(x, p);

        % Complex space vector eccentricity plot (circle or ellipse)
        plot_complex_space_vector(x, t, p);

        % Build scales vector
        [scales, f] = morseScales(p.N, p.fs, p.beta_opt, p.gamma, p.voices_per_octave);

        % Frequency grid and scale spacing
        [da, ~] = scaleSpacing(scales, f);

        % Compute CWT
        fprintf('Running Native 3D Bilateral CWT for [%s]...\n', sig_name);
        [W, dW] = compute_CWT(x, p, scales);

        % Plot bilateral scalogram
        if p.plotting
            plot_CWT_scalogram(W, f, p, t);
        end

        % Apply synchrosqueezing on the CWT
        [Tx_pos, f] = SWT(W(:,:,1), dW(:,:,1), f, scales, p.thresh_ratio, da);
        [Tx_neg, f] = SWT(W(:,:,2), dW(:,:,2), f, scales, p.thresh_ratio, da);

        % Clear intermediate matrix from memory
        clear dW;

        % Normalisation factor for reconstruction and plotting
        C_psi = admissibilityConstant(p.beta_opt, p.gamma);

        % Plot bilateral SWT scalogram
        if p.plot_SWT
            plot_SWT_scalogrm(Tx_pos, Tx_neg, f, C_psi, t, p);
        end

        % Reconstruct signals from CWT and SWT to check integrity
        if p.reconstruct_original_signal

            % Reconstruct signal from bilateral CWT scalogram
            x_hat_cwt = sum(W(:,:,1) .* (da./scales), 1) + sum(W(:,:,2) .* (da./scales), 1);
            x_hat_cwt = x_hat_cwt / C_psi;

            % Reconstruct signal from bilateral SWT scalogram
            x_hat_swt = sum(Tx_pos, 1) + sum(Tx_neg, 1);
            x_hat_swt = x_hat_swt / C_psi;

            plot_reconstructed_signals(t, x, x_hat_cwt, x_hat_swt);

            if p.show_metrics
                metrics_CWT = compute_signal_metrics(x, x_hat_cwt);
                metrics_SWT = compute_signal_metrics(x, x_hat_swt);
            end
        end

        [Tx_all, ~, f_axis, Tx_all_complex] = build_bilateral_matrices(Tx_pos, Tx_neg, W, f, C_psi);

        % Clear redundant memory allocations
        clear W Tx_pos Tx_neg;

        % Energy analysis
        [cyclic_vals2, target_times2, cyclic_vals_rec2, target_times_rec2] = ...
            dual_limit_squeeze_energy(Tx_all, Tx_all_complex, f_axis, t, p, IF, C_psi);

        fprintf('Finished treating [%s]. Freeing allocated memory...\n', sig_name);

        % Clear all signal-specific matrices and variables before next iteration to save RAM
        clear t x currents IF scales f da Tx_all Tx_all_complex ...
              x_hat_cwt x_hat_swt cyclic_vals2 target_times2 cyclic_vals_rec2 target_times_rec2 ...
              metrics_CWT metrics_SWT;
    end
end

function p = params()
    p.fs                                = 5000;
    p.N                                 = 0;
    p.gamma                             = 3;
    p.beta_opt                          = 30;
    p.voices_per_octave                 = 20;
    p.f_disp                            = 50; 
    p.thresh_ratio                      = 1e-4;
    p.plotting                          = true;
    p.plot_SWT                          = true;
    p.plot_rec_sig                      = true;
    p.show_metrics                      = true;
    p.reconstruct_original_signal       = true;
    p.plot_hard_threshold_SWT_direct    = true;
    p.plot_hard_threshold_SWT_rec       = true;
    p.plot_squeeze_SWT_direct           = true;   
    p.plot_squeeze_SWT_rec              = true;   
    p.plot_cyclic_energy                = true;
    p.search_radius_hz                  = 5;
    p.max_expansion_hz                  = 6;
    p.window_sec                        = 0.2;
    
end

%% ==================================================================================================
% PLOT THE THREE PHASE CURRENTS (loaded from the .mat file, before CSV creation)
%===================================================================================================
function plot_three_phase_currents(t, currents, title_str)
    if nargin < 3 || isempty(title_str)
        title_str = 'Three-Phase Currents';
    end
    
    t = t(:);
    ia = currents.ia(:);
    ib = currents.ib(:);
    ic = currents.ic(:);
    
    figure('Color', 'w', 'Name', title_str);
    
    % Overlay plot
    subplot(2,1,1);
    plot(t, ia, 'r', 'LineWidth', 1.0); hold on;
    plot(t, ib, 'g', 'LineWidth', 1.0);
    plot(t, ic, 'b', 'LineWidth', 1.0);
    grid on;
    xlabel('Time (s)'); ylabel('Current (A)');
    title([title_str ' - Overlay']);
    legend('i_a', 'i_b', 'i_c', 'Location', 'best');
    
    % Individual subplots
    ax(1) = subplot(4,1,2);
    plot(t, ia, 'r', 'LineWidth', 1.0); grid on;
    ylabel('i_a (A)');
    title('Phase A');
    
    ax(2) = subplot(4,1,3);
    plot(t, ib, 'g', 'LineWidth', 1.0); grid on;
    ylabel('i_b (A)');
    title('Phase B');
    
    ax(3) = subplot(4,1,4);
    plot(t, ic, 'b', 'LineWidth', 1.0); grid on;
    xlabel('Time (s)');
    ylabel('i_c (A)');
    title('Phase C');
    
    linkaxes(ax, 'x');
end

%% ==================================================================================================
% PLOT DATA FROM A CSV FILE
%===================================================================================================
function plot_csv_data(csv_filename, title_str)
    % Reads a CSV file containing [t, ia, ib, ic] and plots the signals.
    % If the CSV contains a different number of columns, the first column
    % is assumed to be time and all remaining columns are plotted.
    
    if nargin < 2 || isempty(title_str)
        title_str = sprintf('CSV Data: %s', csv_filename);
    end
    
    if ~isfile(csv_filename)
        warning('CSV file not found: %s', csv_filename);
        return;
    end
    
    data = readmatrix(csv_filename);
    
    if isempty(data) || size(data, 2) < 2
        warning('CSV file is empty or has too few columns: %s', csv_filename);
        return;
    end
    
    t = data(:, 1);
    signals = data(:, 2:end);
    num_signals = size(signals, 2);
    
    figure('Color', 'w', 'Name', title_str);
    
    subplot(2, 1, 1);
    hold on; grid on;
    colors = lines(num_signals);
    for k = 1:num_signals
        plot(t, signals(:, k), 'Color', colors(k,:), 'LineWidth', 1.0);
    end
    xlabel('Time (s)'); ylabel('Amplitude');
    title([title_str ' - Overlay']);
    if num_signals == 3
        legend('i_a', 'i_b', 'i_c', 'Location', 'best');
    else
        legend(arrayfun(@(k) sprintf('Col %d', k), 1:num_signals, 'UniformOutput', false), ...
            'Location', 'best');
    end
    
    subplot(2, 1, 2);
    hold on; grid on;
    for k = 1:num_signals
        plot(t, signals(:, k), 'Color', colors(k,:), 'LineWidth', 1.0);
    end
    xlabel('Time (s)'); ylabel('Amplitude');
    title([title_str ' - Zoom (first 0.1 s)']);
    xlim([t(1), t(1) + 0.1]);
end

%% ==================================================================================================
% FFT ANALYSIS AND SPECTRUM PLOT (magnitude only, bilateral, limited to +/-100 Hz)
%===================================================================================================
function plot_FFT_spectrum(x, p)
    % FFT analysis of the complex space vector signal.
    % Plots only the bilateral magnitude spectrum, zoomed to +/-100 Hz.
    
    x = x(:);
    N = length(x);
    
    X = fft(x);
    X = fftshift(X);
    
    % Bilateral frequency axis (works for even and odd N)
    f_axis = (-floor(N/2):ceil(N/2)-1).' * (p.fs / N);
    
    mag = abs(X) / N;
    
    % Limit display to +/-100 Hz
    f_lim = 100;
    keep = (f_axis >= -f_lim) & (f_axis <= f_lim);
    
    figure('Color', 'w', 'Name', 'FFT Magnitude Spectrum');
    plot(f_axis(keep), mag(keep), 'b', 'LineWidth', 1.2);
    grid on;
    xlabel('Frequency (Hz)');
    ylabel('|X(f)|');
    title('FFT Magnitude Spectrum (Bilateral, \pm100 Hz)');
    xlim([-f_lim, f_lim]);
end

%% ==================================================================================================
% COMPLEX SPACE VECTOR ECCENTRICITY (CIRCLE / ELLIPSE)
%===================================================================================================
function plot_complex_space_vector(x, t, p)
    % Complex space vector eccentricity analysis.
    % A healthy machine produces a circular locus; a faulty machine
    % produces an elliptical locus. Eccentricity = sqrt(1 - (b/a)^2).
    
    x = x(:);
    t = t(:);
    
    % Remove DC component
    x_ac = x - mean(x);
    
    % ---- Global ellipse fit via covariance eigen-decomposition ----
    X = [real(x_ac), imag(x_ac)];
    C = cov(X);
    [V, D] = eig(C);
    [eig_vals, sort_idx] = sort(diag(D), 'descend');
    V = V(:, sort_idx);
    
    % Pad to 2 eigenvalues in case the covariance is degenerate (rank 1)
    if numel(eig_vals) < 2
        eig_vals = [eig_vals; 0];
    end
    
    a = 2 * sqrt(max(eig_vals(1), 0));   % semi-major axis
    b = 2 * sqrt(max(eig_vals(2), 0));   % semi-minor axis
    
    if a > 0
        ecc_global = sqrt(1 - (b/a)^2);
    else
        ecc_global = 0;
    end
    
    % ---- Sliding-window eccentricity over time ----
    dt          = t(2) - t(1);
    win_samples = round(p.window_sec / dt);
    N           = length(x);
    num_wins    = ceil(N / win_samples);
    
    ecc_t  = zeros(num_wins, 1);
    t_wins = zeros(num_wins, 1);
    
    for k = 1:num_wins
        i1 = (k-1) * win_samples + 1;
        i2 = min(k * win_samples, N);
        xw = x(i1:i2) - mean(x(i1:i2));
        Xw = [real(xw), imag(xw)];
        Cw = cov(Xw);
        ev = sort(eig(Cw), 'descend');
        
        % Handle degenerate cases (constant data, single point, rank-1 cov)
        if numel(ev) < 2
            ev = [ev; 0];
        end
        a_w = 2 * sqrt(max(ev(1), 0));
        b_w = 2 * sqrt(max(ev(2), 0));
        
        if a_w > 0
            ecc_t(k) = sqrt(max(1 - (b_w/a_w)^2, 0));
        else
            ecc_t(k) = 0;
        end
        t_wins(k) = mean(t(i1:i2));
    end
    
    % ---- Plot ----
    figure('Color', 'w', 'Name', 'Complex Space Vector Eccentricity');
    
    subplot(1, 2, 1);
    plot(real(x_ac), imag(x_ac), 'b.', 'MarkerSize', 2);
    hold on;
    theta = linspace(0, 2*pi, 300);
    if a > 0 && b > 0 && all(isfinite(V(:)))
        ellipse_pts = V * [a*cos(theta); b*sin(theta)];
        plot(ellipse_pts(1,:), ellipse_pts(2,:), 'r-', 'LineWidth', 2);
        legend('Data', 'Fitted Ellipse', 'Location', 'best');
    else
        legend('Data', 'Location', 'best');
    end
    axis equal; grid on;
    xlabel('Real (\alpha)'); ylabel('Imag (\beta)');
    title(sprintf('Complex Space Vector Locus (Ecc = %.4f)', ecc_global));
    
    subplot(1, 2, 2);
    plot(t_wins, ecc_t, 'r-', 'LineWidth', 1.5);
    grid on;
    xlabel('Time (s)'); ylabel('Eccentricity');
    title(sprintf('Eccentricity over Time (%.2fs windows)', p.window_sec));
    ylim([0 1]);
end

%% PLOT BILATERAL SCALOGRAM (CWT)
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
    title("Bilateral Wavelet Scalogram (\pm150 Hz)");
    ylim([-p.f_disp, p.f_disp]); 
    caxis([0 max(amp_all(:))]); colorbar;
    yline(0, '--w', 'LineWidth', 1.5);
    clear amp_all;
end

%% compute_CWT : computes bilateral CWT (positive & negative pages)
function [W, dW] = compute_CWT(x, p, scales)
    x = x(:).';
    Fs = p.fs;
    gamma = p.gamma;
    beta = p.beta_opt;
    
    N = length(x);
    
    X = fft(x);
    w = (2*pi/N) * (0:N-1);
    omega_rad = (2*pi*Fs/N) * [0:floor((N-1)/2), -floor(N/2):-1];
    
    Ns = length(scales);
    
    W  = zeros(Ns, N, 2);
    dW = zeros(Ns, N, 2);
    
    for k = 1:Ns
        a = scales(k);
        % positive page
        psi_pos = morseWaveletFFT(w, a, beta, gamma, 'pos');
        W(k,:,1)  = ifft(X .* conj(psi_pos));
        dW(k,:,1) = ifft(X .* (1i * omega_rad) .* conj(psi_pos));
        % negative page
        psi_neg = morseWaveletFFT(w, a, beta, gamma, 'neg');
        W(k,:,2)  = ifft(X .* conj(psi_neg));
        dW(k,:,2) = ifft(X .* (1i * omega_rad) .* conj(psi_neg));
    end
end

%% morseScales : generates scales and frequency vector
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

%% morseWaveletFFT : frequency-domain Morse wavelet with L1 normalization
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
    
    psi0 = zeros(size(w_domain));
    idx = w_domain > 0;
    
    logA = (beta/gamma) * (1 + log(gamma/beta));
    log_psi0 = logA + beta .* log(w_domain(idx)) - (w_domain(idx)).^gamma;
    psi0(idx) = exp(log_psi0);
    psi = psi0; 
end

%% compute_instantaneous_frequency : estimates IF from CWT coefficients
function IF = compute_instantaneous_frequency(W, dW, thresh_ratio) %#ok<DEFNU>
    if nargin < 3 || isempty(thresh_ratio)
        thresh_ratio = 0.05;
    end
    IF = imag(dW ./ W) / (2*pi);
    thresh = thresh_ratio * max(abs(W), [], 'all');
    IF(abs(W) < thresh) = NaN;
end

%% SWT : synchrosqueezed wavelet transform with integration weight scaling
function [Tx, f_axis] = SWT(W, dW, f, scales, thresh_ratio, da)
    if nargin < 5 || isempty(thresh_ratio), thresh_ratio = 0.05; end
    f = f(:); 
    scales = scales(:);
    [Ns, N] = size(W);

    w_rad = 2 * pi * f;
    [w_sorted, sort_idx] = sort(w_rad, 'ascend');
    row_lookup = (1:Ns)';

    max_abs_W = max(abs(W), [], 'all');
    thresh = thresh_ratio * max_abs_W;

    Tx = zeros(Ns, N, 'like', W);

    for n = 1:N
        W_col = W(:, n);
        dW_col = dW(:, n);

        valid = abs(W_col) > thresh;

        IF_rad = zeros(size(W_col));
        IF_rad(valid) = abs(imag(dW_col(valid) ./ W_col(valid)));

        target_row = nan(Ns, 1);
        if any(valid)
            rows = interp1(w_sorted, row_lookup(sort_idx), IF_rad(valid), 'nearest', 'extrap');
            rows = max(1, min(Ns, rows));
            target_row(valid) = rows;
        end

        valid_indices = find(valid);
        if ~isempty(valid_indices)
            rows = target_row(valid_indices);
            weight_values = W_col(valid_indices) .* (scales(valid_indices) .^ (-1)) .* da(valid_indices);
            Tx_col = accumarray(rows, weight_values, [Ns, 1]);
            Tx(:, n) = Tx_col; 
        end
    end

    f_axis = f;
end

%% scaleSpacing : Delta_a and Delta_omega for a scale/frequency grid
function [da, dw] = scaleSpacing(scales, f)
    scales = scales(:); f = f(:);
    w_rad = 2*pi*f;
    Ns = length(scales);

    da = zeros(Ns,1);
    da(1)       = abs(scales(2)-scales(1));
    da(2:end-1) = abs(scales(3:end)-scales(1:end-2))/2;
    da(end)     = abs(scales(end)-scales(end-1));

    dw = zeros(Ns,1);
    dw(1)       = abs(w_rad(2)-w_rad(1));
    dw(2:end-1) = abs(w_rad(3:end)-w_rad(1:end-2))/2;
    dw(end)     = abs(w_rad(end)-w_rad(end-1));
end

%% admissibilityConstant : closed-form C_psi for the peak-normalized Morse wavelet
function C_psi = admissibilityConstant(beta, gamma) 
    logA = (beta/gamma) * (1 + log(gamma/beta));      
    logC = logA - log(gamma) + gammaln(beta/gamma);
    C_psi = exp(logC);
end

function [t, x] = synthetic_signal(p) 
    Fs = p.fs;                         
    T  = 2;                            
    t  = 0:1/Fs:T-1/Fs;                

    f_pos =  50;                       
    f_neg = -80;                       
    x = exp(1i*2*pi*f_pos*t) + 0.6*exp(1i*2*pi*f_neg*t);
end

function [t, x, currents] = load_faulty_signal(p)
    data = load("faulty_stator_v00.mat"); 
    N = length(data.Te_f_12);
    fs=p.fs;
    t = (0:N-1) / fs;
    idx = (t >= 2);
    t = t(idx);

    ia = data.isa_f_12(idx, 1);
    ib = data.isb_f_12(idx, 1);
    ic = data.isc_f_12(idx, 1);

    x = sqrt(2/3) * (ia + ib * exp(1i*2*pi/3) + ic * exp(1i*4*pi/3));
    
    % Return raw three-phase currents for standalone plotting
    currents = struct('ia', ia, 'ib', ib, 'ic', ic);
end

function [t, x, currents] = load_healthy_signal(p)
    data = load("healthy_v00.mat"); 
    N = length(data.Te_h);
    fs=p.fs;
    t = (0:N-1) / fs;
    idx = (t >= 2);
    t = t(idx);

    ia = data.isa_h(idx, 1);
    ib = data.isb_h(idx, 1);
    ic = data.isc_h(idx, 1);

    x = sqrt(2/3) * (ia + ib * exp(1i*2*pi/3) + ic * exp(1i*4*pi/3));
    
    % Return raw three-phase currents for standalone plotting
    currents = struct('ia', ia, 'ib', ib, 'ic', ic);
end

function [t, x] = load_chirp_signal(p) 
    fs = p.fs;              
    t  = (0:1/fs:2-1/fs)';  

    f0 = 50; f1 = 200; T = t(end);
    k  = (f1 - f0) / T;     
    phase = 2*pi*(f0*t + k/2*t.^2);
    x1 = exp(1i*phase);     

    x2 = 0.6*exp(-1i*2*pi*80*t); 

    x = x1 + x2;            
    x = x + 0.02*(randn(size(x)) + 1i*randn(size(x)));
end

function plot_SWT_scalogrm(Tx_pos, Tx_neg, f, C_psi, t, p) 
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
    clear amp_all_swt;
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

    p = abs(x).^2;
    p = p / sum(p);
    p(p==0)=[];

    metrics.Shannon_entropy = -sum(p .* log2(p));

    alpha = 2;
    metrics.Renyi_entropy = 1/(1-alpha) * log2(sum(p.^alpha));

    x_real = real(x);
    xhat_real = real(x_hat);

    C = corrcoef(x_real,xhat_real);
    metrics.Correlation = C(1,2);

    fprintf('\n========== Signal Metrics ==========\n')
    fprintf('RMSE                  : %.6e\n', metrics.RMSE);
    fprintf('NRMSE                 : %.6e\n', metrics.NRMSE);
    fprintf('SNR                   : %.3f dB\n', metrics.SNR);
    fprintf('Energy original       : %.6e\n', metrics.Energy_original);
    fprintf('Energy reconstructed  : %.6e\n', metrics.Energy_reconstructed);
    fprintf('Shannon entropy       : %.6f bits\n', metrics.Shannon_entropy);
    fprintf('Renyi entropy         : %.6f bits\n', metrics.Renyi_entropy);
    fprintf('Correlation = %.12f\n',metrics.Correlation)
    fprintf('====================================\n')
end

function [t, x, IF] = load_signal_freq()
    fs        = 5000;      
    T_total   = 4;          
    f0        = 20;         
    f1        = 30;         
    T_const   = 1;          
    T_trans   = T_total - 2*T_const;  
    A         = 0.3;        

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

%% ==================================================================================================
% First threshold METHOD: Hard threshold
%===================================================================================================
function [cyclic_vals1, target_times1, cyclic_vals_rec1, target_times_rec1] = hard_threshold_energy(Tx_all, Tx_all_complex, W_all, f_axis, t, p, IF, C_psi)

    % 1. Estimate global noise threshold using Donoho's Universal Threshold
    t_thresh = estimate_donoho_threshold(W_all, length(t))
    
    % 2. Extract adaptive frequency bands and compute instantaneous energy
    [band_min_f, ~, band_max_f_vis, instant_energy, band_indices] = track_adaptive_bands(...
        Tx_all, f_axis, IF, t_thresh, p.search_radius_hz, p.max_expansion_hz);

    % Method 1: compute directly from the SWT scalogram
    [cyclic_time, cyclic_energy] = compute_cyclic_energy(instant_energy, t, p.window_sec);
    
    if p.plot_hard_threshold_SWT_direct
        plot_ridge_energy_results(t, f_axis, IF, Tx_all, band_min_f, band_max_f_vis, ...
            instant_energy, cyclic_time, cyclic_energy, p);
    end
    
    [cyclic_vals1, target_times1] = track_cyclic_energy(instant_energy, t, p);
    
    % Method 2: compute energy from signal reconstruction
    [x_rec, instant_energy_rec] = reconstruct_signal_from_band(...
        Tx_all_complex, Tx_all, band_indices, t_thresh, C_psi);
    
    clear Tx_all_complex; % Free complex spectrum matrix

    [cyclic_time, cyclic_energy_rec] = compute_cyclic_energy(instant_energy_rec, t, p.window_sec);

    if p.plot_hard_threshold_SWT_rec
        plot_reconstruction_results(t, f_axis, IF, Tx_all, band_min_f, band_max_f_vis, ...
            x_rec, instant_energy_rec, cyclic_time, cyclic_energy_rec, p);
    end
    
    [cyclic_vals_rec1, target_times_rec1] = track_cyclic_energy(instant_energy_rec, t, p);
end

function [Tx_all, W_all, f_axis, Tx_all_complex] = build_bilateral_matrices(Tx_pos, Tx_neg, W, f, C_psi)
    NumTime = size(Tx_pos, 2);
    
    f_axis = [-f; 0; flipud(f)];
    [f_axis, sortIdx] = sort(f_axis, 'ascend');

    % Magnitudes (for visualization & thresholding)
    Tx_all = [abs(Tx_neg); zeros(1, NumTime); flipud(abs(Tx_pos))] / C_psi;
    Tx_all = Tx_all(sortIdx, :);

    % Complex matrix (for accurate reconstruction)
    Tx_all_complex = [Tx_neg; zeros(1, NumTime); flipud(Tx_pos)];
    Tx_all_complex = Tx_all_complex(sortIdx, :);

    W_all = [abs(W(:,:,2)); zeros(1, NumTime); flipud(abs(W(:,:,1)))] / C_psi;
    W_all = W_all(sortIdx, :);
end

function t_thresh = estimate_donoho_threshold(W_all, N)
    highest_freq_row = W_all(end, :);
    mad_val = median(abs(highest_freq_row - median(highest_freq_row)));  
    sigma_hat = mad_val / 0.6745;
    %t_thresh = (sigma_hat / sqrt(N)) * sqrt(2 * log(N - 1));
    t_thresh = sigma_hat * sqrt(2*log(N));
end

function [band_min_f, band_max_f, band_max_f_vis, instant_energy, band_indices] = track_adaptive_bands(...
    Tx_all, f_axis, IF, t_thresh, search_radius_hz, max_expansion_hz)
    
    [NumFreq, NumTime] = size(Tx_all);
    
    band_min_f     = zeros(NumTime, 1);
    band_max_f     = zeros(NumTime, 1);
    band_max_f_vis = zeros(NumTime, 1);
    instant_energy = zeros(NumTime, 1);
    band_indices   = cell(NumTime, 1);

    for idx = 1:NumTime
        current_IF = IF(idx);
        [~, nominal_f_idx] = min(abs(f_axis - current_IF));

        search_range = find(f_axis >= (current_IF - search_radius_hz) & ...
                            f_axis <= (current_IF + search_radius_hz));

        if isempty(search_range)
            search_range = nominal_f_idx;
        end
        
        [~, local_peak_offset] = max(Tx_all(search_range, idx));
        center_f_idx = search_range(local_peak_offset);
        f_center     = f_axis(center_f_idx);

        limit_up   = find(f_axis <= (f_center + max_expansion_hz), 1, 'last');
        limit_down = find(f_axis >= (f_center - max_expansion_hz), 1, 'first');

        if isempty(limit_up),   limit_up   = NumFreq; end
        if isempty(limit_down), limit_down = 1;       end

        up_idx = center_f_idx;
        while up_idx < limit_up && Tx_all(up_idx + 1, idx) >= t_thresh
            up_idx = up_idx + 1;
        end

        down_idx = center_f_idx;
        while down_idx > limit_down && Tx_all(down_idx - 1, idx) >= t_thresh
            down_idx = down_idx - 1;
        end

        down_idx = max(1, down_idx);
        up_idx   = min(NumFreq, up_idx);

        band_min_f(idx) = f_axis(down_idx);
        band_max_f(idx) = f_axis(up_idx);
        
        if up_idx < NumFreq
            band_max_f_vis(idx) = f_axis(up_idx + 1);
        else
            band_max_f_vis(idx) = f_axis(up_idx);
        end

        band_range        = down_idx:up_idx;
        band_indices{idx} = band_range;

        valid_coefs = Tx_all(band_range, idx);
        denoised_coefs = valid_coefs;
        denoised_coefs(valid_coefs < t_thresh) = 0;

        instant_energy(idx) = sum(denoised_coefs.^2);
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

function plot_ridge_energy_results(t, f_axis, IF, Tx_all, band_min_f, band_max_f_vis, ...
    instant_energy, cyclic_time, cyclic_energy, p)
    
    keep = (f_axis >= -p.f_disp) & (f_axis <= p.f_disp);
    
    figure('Color', 'w', 'Name', 'SWT Scalogram & Ridge Tracking');
    h_surf = pcolor(t, f_axis(keep), Tx_all(keep, :));
    set(h_surf, 'EdgeColor', 'none', 'FaceColor', 'flat');
    axis xy; hold on;
    
    z_top = ones(size(t)) * max(Tx_all(:)) * 2;
    h_if   = plot3(t, IF, z_top, 'w--', 'LineWidth', 2);
    h_band = plot3(t, band_min_f, z_top, 'r-', 'LineWidth', 1.5);
    plot3(t, band_max_f_vis, z_top, 'r-', 'LineWidth', 1.5);
    
    xlabel('Time (s)'); ylabel('Frequency (Hz)');
    title('SWT Scalogram with Selected Adaptive Band & True IF');
    ylim([-p.f_disp, p.f_disp]); colorbar; caxis([0 max(Tx_all(:))]);
    legend([h_if, h_band], {'True IF Profile', 'Adaptive Band Boundaries'}, 'Location', 'Best');
    hold off;

    figure('Color', 'w', 'Name', 'Energy Profiles');
    subplot(2,1,1); 
    plot(t, instant_energy, 'b', 'LineWidth', 1.5); grid on;
    xlabel('Time (s)'); ylabel('Energy (Pure SWT domain)');
    title('Instantaneous Energy inside Tracked SWT Band');
    
    subplot(2,1,2); 
    stem(cyclic_time, cyclic_energy, 'r', 'Filled', 'LineWidth', 1.5); grid on;xlim([2 4]);
    xlabel('Time (s)'); ylabel('Accumulated Energy'); 
    title('Cyclic Energy Profile (Every 0.2 seconds)');
end

function [x_rec, instant_energy_rec] = reconstruct_signal_from_band(...
    Tx_all_complex, Tx_all, band_indices, t_thresh, C_psi)
    
    NumTime = length(band_indices);
    x_rec = zeros(NumTime, 1);
    instant_energy_rec = zeros(NumTime, 1);

    for idx = 1:NumTime
        band_range = band_indices{idx};
        
        coefs_tx_complex = Tx_all_complex(band_range, idx);
        valid_swt_check  = Tx_all(band_range, idx);

        coefs_tx_complex(valid_swt_check < t_thresh) = 0;

        x_rec(idx) = sum(coefs_tx_complex, 1) / C_psi;

        instant_energy_rec(idx) = abs(x_rec(idx))^2;
    end
end

function plot_reconstruction_results(t, f_axis, IF, Tx_all, band_min_f, band_max_f_vis, ...
    x_rec, instant_energy_rec, cyclic_time, cyclic_energy_rec, p)
    
    keep = (f_axis >= -p.f_disp) & (f_axis <= p.f_disp);

    figure('Color', 'w', 'Name', 'SWT Scalogram & Band Tracking');
    h_surf = pcolor(t, f_axis(keep), Tx_all(keep, :));
    set(h_surf, 'EdgeColor', 'none', 'FaceColor', 'flat');
    axis xy; hold on;

    z_top = ones(size(t)) * max(Tx_all(:)) * 2;
    h_if   = plot3(t, IF, z_top, 'w--', 'LineWidth', 2);
    h_band = plot3(t, band_min_f, z_top, 'r-', 'LineWidth', 1.5);
    plot3(t, band_max_f_vis, z_top, 'r-', 'LineWidth', 1.5);

    xlabel('Time (s)'); ylabel('Frequency (Hz)');
    title('SWT Scalogram with Selected Adaptive Band & SWT Inversion');
    ylim([-p.f_disp, p.f_disp]); colorbar; caxis([0 max(Tx_all(:))]);
    legend([h_if, h_band], {'True IF Profile', 'Adaptive Band Boundaries'}, 'Location', 'Best');
    hold off;

    figure('Color', 'w', 'Name', 'Reconstructed Signal & Energy Profiles');
    subplot(3,1,1); 
    plot(t, real(x_rec), 'g', 'LineWidth', 1.2); grid on;
    xlabel('Time (s)'); ylabel('Amplitude'); 
    title('Denoised Reconstructed Signal (Real Part)');

    subplot(3,1,2); 
    plot(t, instant_energy_rec, 'b', 'LineWidth', 1.5); grid on;
    xlabel('Time (s)'); ylabel('Instantaneous Energy'); 
    title('Instantaneous Temporal Energy inside Tracked Band');

    subplot(3,1,3); 
    stem(cyclic_time, cyclic_energy_rec, 'r', 'Filled', 'LineWidth', 1.5); grid on;
    xlabel('Time (s)'); ylabel('Accumulated Energy'); 
    title('Cyclic Energy Profile of Reconstructed Signal (Every 0.2 seconds)');
end

%% ==================================================================================================
% Second threshold METHOD: Dual-Limit Squeeze Threshold
%===================================================================================================
function [cyclic_vals2, target_times2, cyclic_vals_rec2, target_times_rec2] = dual_limit_squeeze_energy(Tx_all, Tx_all_complex, f_axis, t, p, IF, C_psi)

    % 1. Extract adaptive frequency bands and compute instantaneous energy using dual-limit squeeze
    [band_min_f, ~, band_max_f_vis, instant_energy, band_indices, uth_history, lth_history, ...
        iteration_counts, sample_indices, convergence_paths] = track_adaptive_squeeze_bands(...
        Tx_all, f_axis, IF, p.search_radius_hz, p.max_expansion_hz);

    % Method 1: compute directly from the SWT scalogram
    [cyclic_time, cyclic_energy] = compute_cyclic_energy(instant_energy, t, p.window_sec);

    if p.plot_squeeze_SWT_direct
        plot_squeeze_direct_results(t, f_axis, IF, Tx_all, band_min_f, band_max_f_vis, ...
            instant_energy, cyclic_time, cyclic_energy, uth_history, lth_history, ...
            iteration_counts, sample_indices, convergence_paths, p);
    end
    
    [cyclic_vals2, target_times2] = track_cyclic_energy(instant_energy, t, p);
    

    % Method 2: compute energy from signal reconstruction
    [x_rec, instant_energy_rec] = reconstruct_squeeze_signal_from_band(...
        Tx_all_complex, Tx_all, band_indices, lth_history, C_psi);

    clear Tx_all_complex; % Free complex spectrum matrix

    [cyclic_time, cyclic_energy_rec] = compute_cyclic_energy(instant_energy_rec, t, p.window_sec);

    if p.plot_squeeze_SWT_rec
        plot_squeeze_reconstruction_results(t, f_axis, IF, Tx_all, band_min_f, band_max_f_vis, ...
            x_rec, instant_energy_rec, cyclic_time, cyclic_energy_rec, uth_history, lth_history, ...
            sample_indices, convergence_paths, p);
    end
    
    [cyclic_vals_rec2, target_times_rec2] = track_cyclic_energy(instant_energy_rec, t, p);
    
end

function [band_min_f, band_max_f, band_max_f_vis, instant_energy, band_indices, ...
    uth_history, lth_history, iteration_counts, sample_indices, convergence_paths] = track_adaptive_squeeze_bands(...
    Tx_all, f_axis, IF, search_radius_hz, max_expansion_hz)

    [NumFreq, NumTime] = size(Tx_all);

    band_min_f      = zeros(NumTime, 1);
    band_max_f      = zeros(NumTime, 1);
    band_max_f_vis  = zeros(NumTime, 1);
    instant_energy  = zeros(NumTime, 1);
    band_indices    = cell(NumTime, 1);
    uth_history     = zeros(NumTime, 1);
    lth_history     = zeros(NumTime, 1);

    w = 0.125;           
    max_iters = 100;     

    iteration_counts = zeros(NumTime, 1);
    sample_indices = [round(0.10 * NumTime), round(0.50 * NumTime), round(0.90 * NumTime)];
    convergence_paths = struct();

    for idx = 1:NumTime
        current_IF = IF(idx);
        [~, nominal_f_idx] = min(abs(f_axis - current_IF));

        search_range = find(f_axis >= (current_IF - search_radius_hz) & ...
                            f_axis <= (current_IF + search_radius_hz));

        if isempty(search_range)
            search_range = nominal_f_idx;
        end
        slice_amps = Tx_all(search_range, idx);

        [a, local_peak_offset] = max(slice_amps);
        center_f_idx = search_range(local_peak_offset);
        f_center     = f_axis(center_f_idx);

        if a < eps
            band_min_f(idx)     = f_center;
            band_max_f(idx)     = f_center;
            band_max_f_vis(idx) = f_axis(min(NumFreq, center_f_idx + 1));
            band_indices{idx}   = center_f_idx;
            continue;
        end

        uth = 0.50 * a;             %upper threshold initialization
        lth = 0.005 * a;             % lower threshold inintialization

        is_sampled = ismember(idx, sample_indices);
        if is_sampled
            sample_tag = sprintf('idx_%d', idx);
            convergence_paths.(sample_tag).uth = uth;
            convergence_paths.(sample_tag).lth = lth;
        end

        actual_iters = max_iters;
        for iter = 1:max_iters
            D = uth - lth;
            uth = uth - (w * D);
            lth = lth + (w * D);

            N_uth = count_peaks_above_thresh(slice_amps, uth);
            N_lth = count_peaks_above_thresh(slice_amps, lth);

            if is_sampled
                convergence_paths.(sample_tag).uth(end+1) = uth;
                convergence_paths.(sample_tag).lth(end+1) = lth;
            end

            if N_uth == N_lth || (uth - lth) < 1e-6 * a
                actual_iters = iter;
                
                break;
            end
        end
        iteration_counts(idx) = actual_iters;

        uth_history(idx) = uth;
        lth_history(idx) = lth;
        
        limit_up   = find(f_axis <= (current_IF + max_expansion_hz), 1, 'last');
        limit_down = find(f_axis >= (current_IF - max_expansion_hz), 1, 'first');

        if isempty(limit_up),   limit_up   = NumFreq; end
        if isempty(limit_down), limit_down = 1;       end

        up_idx = center_f_idx;
        while up_idx < limit_up && Tx_all(up_idx + 1, idx) >= lth
            up_idx = up_idx + 1;
        end

        down_idx = center_f_idx;
        while down_idx > limit_down && Tx_all(down_idx - 1, idx) >= lth
            down_idx = down_idx - 1;
        end

        down_idx = max(1, down_idx);
        up_idx   = min(NumFreq, up_idx);

        band_min_f(idx) = f_axis(down_idx);
        band_max_f(idx) = f_axis(up_idx);

        if up_idx < NumFreq
            band_max_f_vis(idx) = f_axis(up_idx + 1);
        else
            band_max_f_vis(idx) = f_axis(up_idx);
        end

        band_range        = down_idx:up_idx;
        band_indices{idx} = band_range;

        valid_coefs = Tx_all(band_range, idx);
        denoised_coefs = valid_coefs;
        denoised_coefs(valid_coefs < lth) = 0;

        instant_energy(idx) = sum(denoised_coefs.^2);
    end
end

function num_peaks = count_peaks_above_thresh(signal_slice, thresh)
    above_idx = find(signal_slice > thresh);
    if isempty(above_idx)
        num_peaks = 0;
        return;
    end

    is_peak = false(size(signal_slice));
    for i = 2:(length(signal_slice) - 1)
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
    Tx_all_complex, Tx_all, band_indices, lth_history, C_psi)

    NumTime = length(band_indices);
    x_rec = zeros(NumTime, 1);
    instant_energy_rec = zeros(NumTime, 1);

    for idx = 1:NumTime
        band_range = band_indices{idx};
        lth = lth_history(idx);

        coefs_tx_complex = Tx_all_complex(band_range, idx);
        valid_swt_check  = Tx_all(band_range, idx);

        coefs_tx_complex(valid_swt_check < lth) = 0;

        x_rec(idx) = sum(coefs_tx_complex, 1) / C_psi;
        instant_energy_rec(idx) = abs(x_rec(idx))^2;
    end
end

function plot_squeeze_direct_results(t, f_axis, IF, Tx_all, band_min_f, band_max_f_vis, ...
    instant_energy, cyclic_time, cyclic_energy, uth_history, lth_history, ...
    iteration_counts, sample_indices, convergence_paths, p)

    keep = (f_axis >= -p.f_disp) & (f_axis <= p.f_disp);

    % Figure 1: Scalogram Overlay
    figure('Color', 'w', 'Name', 'Adaptive Squeeze Scalogram Tracking');
    h_surf = pcolor(t, f_axis(keep), Tx_all(keep, :));
    set(h_surf, 'EdgeColor', 'none', 'FaceColor', 'flat');
    axis xy; hold on;

    z_top = ones(size(t)) * max(Tx_all(:)) * 2;
    h_if   = plot3(t, IF, z_top, 'w--', 'LineWidth', 2);
    h_band = plot3(t, band_min_f, z_top, 'r-', 'LineWidth', 1.5);
    plot3(t, band_max_f_vis, z_top, 'r-', 'LineWidth', 1.5);

    xlabel('Time (s)'); ylabel('Frequency (Hz)');
    title('SWT Scalogram with Adaptive Squeeze Band Tracking');
    ylim([-p.f_disp, p.f_disp]); colorbar; caxis([0 max(Tx_all(:))*0.01]);
    legend([h_if, h_band], {'True IF Profile', 'Adaptive Squeeze Boundaries'}, 'Location', 'Best');
    hold off;

    % Figure 2: Energy Profiles
    figure('Color', 'w', 'Name', 'Energy Profiles');
    subplot(2,1,1);
    plot(t, instant_energy, 'b', 'LineWidth', 1.5); grid on;
    xlabel('Time (s)'); ylabel('Energy (Pure SWT domain)');
    title('Instantaneous Energy inside Tracked Band');

    subplot(2,1,2);
    stem(cyclic_time, cyclic_energy, 'r', 'Filled', 'LineWidth', 1.5); grid on;
    xlabel('Time (s)'); ylabel('Accumulated Energy');
    title('Cyclic Energy Profile (Every 0.2 seconds)');

   
end

function plot_squeeze_reconstruction_results(t, f_axis, IF, Tx_all, band_min_f, band_max_f_vis, ...
    x_rec, instant_energy_rec, cyclic_time, cyclic_energy_rec, uth_history, lth_history, ...
    sample_indices, convergence_paths, p)

    keep = (f_axis >= -p.f_disp) & (f_axis <= p.f_disp);

    % Figure 1: Scalogram Overlay
    figure('Color', 'w', 'Name', 'SWT Scalogram & Squeeze Band');
    h_surf = pcolor(t, f_axis(keep), Tx_all(keep, :));
    set(h_surf, 'EdgeColor', 'none', 'FaceColor', 'flat');
    axis xy; hold on;

    z_top = ones(size(t)) * max(Tx_all(:)) * 2;
    h_if   = plot3(t, IF, z_top, 'w--', 'LineWidth', 2);
    h_band = plot3(t, band_min_f, z_top, 'r-', 'LineWidth', 1.5);
    plot3(t, band_max_f_vis, z_top, 'r-', 'LineWidth', 1.5);

    xlabel('Time (s)'); ylabel('Frequency (Hz)');
    title('SWT Scalogram with Selected Adaptive Band & SWT Inversion');
    ylim([-p.f_disp, p.f_disp]); colorbar; caxis([0 max(Tx_all(:))*0.01]);
    legend([h_if, h_band], {'True IF Profile', 'Adaptive Band Boundaries'}, 'Location', 'Best');
    hold off;

    % Figure 2: Signal Reconstruction & Energy Profiles
    figure('Color', 'w', 'Name', 'Reconstructed Signal & Temporal Energy');
    subplot(3,1,1);
    plot(t, real(x_rec), 'g', 'LineWidth', 1.2); grid on;
    xlabel('Time (s)'); ylabel('Amplitude');
    title('Denoised Reconstructed Signal (Real Part)'); xlim([0.5 1]);

    subplot(3,1,2);
    plot(t, instant_energy_rec, 'b', 'LineWidth', 1.5); grid on;
    xlabel('Time (s)'); ylabel('Instantaneous Energy');
    title('Instantaneous Temporal Energy inside Tracked Band');

    subplot(3,1,3);
    stem(cyclic_time, cyclic_energy_rec, 'r', 'Filled', 'LineWidth', 1.5); grid on;
    xlabel('Time (s)'); ylabel('Accumulated Energy');
    title('Cyclic Energy Profile of Reconstructed Signal (Every 0.2 seconds)');

    
end

function signal_noisy = add_custom_wgn(signal_clean, target_snr_db)
    % 1. Compute the true power (variance) of your clean signal
    % (Works perfectly for any signal shape of amplitude 0.3)
    signal_power = mean(signal_clean.^2);
    
    % 2. Convert target SNR from dB to linear scale
    snr_linear = 10^(target_snr_db / 10);
    
    % 3. Calculate the required noise power (variance)
    noise_power = signal_power / snr_linear;
    
    % 4. Calculate standard deviation (rms amplitude of the noise)
    sigma = sqrt(noise_power);
    
    % 5. Generate Zero-Mean White Gaussian Noise and add it
    noise = sigma * randn(size(signal_clean));
    signal_noisy = signal_clean + noise;
end

function [cyclic_vals, target_times] = track_cyclic_energy(instant_energy, t, p, method_name)
    % 1. Determine discrete time targets starting from the minimum time vector limit
    t_start = t(1);
    t_end   = t(end);
    target_times = (t_start : p.window_sec : t_end).';
    
    % 2. Use nearest neighbor matching to pull exact instantaneous values
    num_points = length(target_times);
    cyclic_vals = zeros(num_points, 1);
    
    for i = 1:num_points
        [~, exact_idx] = min(abs(t - target_times(i)));
        cyclic_vals(i) = instant_energy(exact_idx);
    end
    
    % 3. Check flag parameter configuration to dynamically plot
    if p.plot_cyclic_energy
        figure;
        plot(target_times, cyclic_vals, 'LineWidth', 1.5); 
        grid on;
        xlabel('Time (s)'); 
        ylabel('Instantaneous Energy Value');
        title("Cyclic Energy Tracking");
    end
end