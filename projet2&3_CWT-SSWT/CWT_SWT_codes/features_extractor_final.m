function main()

    t_total = tic;

    p = config();

    script_dir = fileparts(mfilename('fullpath'));
    if isempty(script_dir), script_dir = pwd; end

    sample1 = load_healthy_sample(p);
    sample1.label = 'Healthy';

    sample2 = load_faulty_sample(p);
    sample2.label = 'Faulty';

    db = [sample1, sample2];

    n_samples = length(db);
    feature_names = {};
    sample_feature_evolution = cell(1, n_samples); % {k} -> [n_windows_k x n_features]
    window_times_all         = cell(1, n_samples); % {k} -> [n_windows_k x 1]

    for k = 1:n_samples
        s = db(k);

        t_step = tic; % Start sample timer before extraction

        % Cyclic windowed extraction over the whole recording
        [feat_evo, f_names, win_times] = extract_features_cyclic(s, p);

        t_elapsed = toc(t_step); % Stop sample timer right after extraction

        sample_feature_evolution{k} = feat_evo;
        window_times_all{k}         = win_times;
        feature_names                = f_names;

        fprintf("  -> Sample %d (%s) took : %.3f seconds [%d windows of %.3f s]\n", ...
            k, s.label, t_elapsed, size(feat_evo, 1), p.window_length);
    end

    % Reorganize results: one struct per feature, holding its evolution per sample
    feature_struct_array = build_feature_structs(feature_names, sample_feature_evolution, window_times_all, db);

    % Relative variation (%) of each feature between Healthy and Faulty,
    % sorted by discriminative power (largest |variation| first)
    relvar_table = compute_relative_variation(feature_struct_array);
    
    % Display table
    disp(relvar_table);

    % Save CSV table
    relvar_filepath = fullfile(script_dir, 'relative_variance.csv');
    writetable(relvar_table, relvar_filepath);
    fprintf('Relative variation table saved to:\n%s\n', relvar_filepath);

    % 1. Plot top 25 features for mentor presentation slides (clean and highly legible)
    plot_top_relative_variance_barchart(relvar_table, 50);

    % 2. Plot all 117 features in a tall reference bar chart
    plot_relative_variance_barchart(relvar_table);

%     % Plot the evolution of every feature, one figure per feature (on screen)
%     plot_feature_evolutions(feature_struct_array);
% 
%     % Export the very same figures into a single multi-page PDF report
%     report_filepath = fullfile(script_dir, 'feature_evolution_report.pdf');
%     export_figures_to_pdf(feature_struct_array, report_filepath);
% 
%     % Save candidate feature pool (windowed) to the script's directory
%     output_filepath = fullfile(script_dir, 'candidate_feature_pool_windowed.mat');
% 
%     save(output_filepath, 'feature_struct_array', 'feature_names', 'relvar_table', 'db', 'p');
%     fprintf('Extraction complete: %d features x %d samples saved to:\n%s\n', ...
%             length(feature_names), n_samples, output_filepath);
%     fprintf("\nDone! Total execution time: %.3f seconds\n", toc(t_total));
end

%% CONFIGURATION
function p = config()
    p.fs = 5000;            % Sampling frequency (Hz)
    p.f0 = 50;              % Fundamental grid frequency (Hz)
    p.H  = 5;               % Maximum harmonic order
    p.window_length = 0.2;  % Analysis window length (s) — default
end

%% WINDOW SEGMENTATION (cyclic slicing of an already-computed continuous sample)
function windows = segment_into_windows(s, p)
    window_samples = round(p.window_length * p.fs);
    N = length(s.i_pos);
    n_windows = floor(N / window_samples);

    windows = repmat(struct('t', [], 'is', [], 'i_zero', [], 'i_pos', [], ...
                            'i_neg', [], 'Te', [], 't_start', 0), 1, n_windows);

    for w = 1:n_windows
        idx_range = (w-1)*window_samples + 1 : w*window_samples;

        windows(w).t       = s.t(idx_range);
        windows(w).is      = s.is(idx_range);
        windows(w).i_zero  = s.i_zero(idx_range);
        windows(w).i_pos   = s.i_pos(idx_range);
        windows(w).i_neg   = s.i_neg(idx_range);
        windows(w).Te      = s.Te(idx_range);
        windows(w).t_start = s.t(idx_range(1));
    end
end

%% CYCLIC FEATURE EXTRACTION OVER WINDOWS
function [feat_evo, f_names, win_times] = extract_features_cyclic(s, p)
    windows   = segment_into_windows(s, p);
    n_windows = length(windows);

    feat_evo  = [];
    f_names   = {};
    win_times = zeros(n_windows, 1);

    for w = 1:n_windows
        [f_vec, names] = extract_all_features(windows(w), p);
        feat_evo = [feat_evo; f_vec]; %#ok<AGROW>
        f_names  = names;
        win_times(w) = windows(w).t_start;
    end
end

%% FEATURE EXTRACTION MASTER PIPELINE
function [f_vec, f_names] = extract_all_features(s, p)
    % 1. Positive Sequence Component (I1)
    [I1_time, time_labels] = extract_time_features(s.i_pos);
    [I1_freq, freq_labels] = extract_freq_features(s.i_pos, p.fs, p.f0, p.H);
    
    % 2. Negative Sequence Component (I2)
    [I2_time, ~]           = extract_time_features(s.i_neg);
    [I2_freq, ~]           = extract_freq_features(s.i_neg, p.fs, p.f0, p.H);
    
    % 3. Cross-Component (Derived from I1, I2)
    cuf_val                = compute_cuf(s.i_pos, s.i_neg);
    
    % 4. Complex Space Vector (CSV / i_qds)
    [CSV_time, ~]          = extract_time_features(s.is);
    [CSV_freq, ~]          = extract_freq_features(s.is, p.fs, p.f0, p.H);
    [CSV_geom, geom_labels]= extract_geometric_features(s.is);
    [CSV_nonl, nonl_labels]= extract_nonlinear_features(s.is);
    
    % 5. Electromagnetic Torque (Te)
    [Te_time, ~]           = extract_time_features(s.Te);
    
    % Concatenate full feature vector
    f_vec = [ ...
        I1_time, I1_freq, ...
        I2_time, I2_freq, ...
        cuf_val, ...
        CSV_time, CSV_freq, CSV_geom, CSV_nonl, ...
        Te_time ...
    ];

    % Generate prefix-scoped feature column names
    f_names = [ ...
        strcat('I1_', time_labels), strcat('I1_', freq_labels), ...
        strcat('I2_', time_labels), strcat('I2_', freq_labels), ...
        {'Cross_CUF'}, ...
        strcat('CSV_', time_labels), strcat('CSV_', freq_labels), ...
        strcat('CSV_', geom_labels), strcat('CSV_', nonl_labels), ...
        strcat('Te_', time_labels) ...
    ];
end

%% DOMAIN EXTRACTION HELPERS
function [f_time, labels] = extract_time_features(x)
    f_time = [ ...
        compute_peak(x), compute_rms(x), compute_kurtosis(x), ...
        compute_crest_factor(x), compute_clearance_factor(x), compute_impulse_factor(x), ...
        compute_shape_factor(x), compute_skewness(x), compute_square_mean_root(x), ...
        compute_5th_norm_moment(x), compute_6th_norm_moment(x), compute_mean(x), ...
        compute_shape_factor_2(x), compute_peak_to_peak(x), compute_kurtosis_factor(x), ...
        compute_std(x), compute_smoothness(x), compute_uniformity(x), ...
        compute_neg_log_likelihood(x), compute_ac_dc_ratio(x), compute_tkeo_mean(x), ...
        compute_tkeo_kurtosis(x) ...
    ];
    labels = { ...
        'Peak', 'RMS', 'Kurtosis', 'Crest_Factor', 'Clearance_Factor', ...
        'Impulse_Factor', 'Shape_Factor', 'Skewness', 'SMR', 'M5_Norm', ...
        'M6_Norm', 'Mean', 'Shape_Factor_2', 'Peak_to_Peak', 'Kurtosis_Factor', ...
        'Std', 'Smoothness', 'Uniformity', 'Neg_Log_Likelihood', 'AC_DC_Ratio', ...
        'TKEO_Mean', 'TKEO_Kurtosis' ...
    };
end

function [f_freq, labels] = extract_freq_features(x, fs, f0, H)
    f_freq = [ ...
        compute_spectral_centroid(x, fs), ...
        compute_spectral_spread(x, fs), ...
        compute_spectral_kurtosis(x, fs), ...
        compute_target_band_energy(x, fs, 2*f0-10, 2*f0+10), ...
        compute_harmonic_distortion(x, fs, f0, H), ...
        compute_spectral_flatness(x, fs), ...
        compute_kurtogram_peak(x, fs) ...
    ];
    labels = { ...
        'Spectral_Centroid', 'Spectral_Spread', 'Spectral_Kurtosis', ...
        'Target_Band_Energy', 'Harmonic_Distortion', 'Spectral_Flatness', 'Kurtogram_Peak' ...
    };
end

function [f_geom, labels] = extract_geometric_features(is)
    f_geom = [ ...
        compute_eccentricity(is), ...
        compute_orbit_thickness(is), ...
        compute_enclosed_area(is) ...
    ];
    labels = {'Eccentricity', 'Orbit_Thickness', 'Enclosed_Area'};
end

function [f_nonl, labels] = extract_nonlinear_features(x)
    f_nonl = [ ...
        compute_sample_entropy(x, 2, 0.2*std(abs(x))), ...
        compute_higuchi_fractal_dim(x, 10), ...
        compute_lempel_ziv_complexity(x), ...
        compute_hilbert_envelope_entropy(x) ...
    ];
    labels = {'Sample_Entropy', 'Higuchi_Fractal_Dim', 'Lempel_Ziv_Complexity', 'Hilbert_Envelope_Entropy'};
end

%% DATA LOADING
function signal = load_faulty_sample(p)
    data = load("faulty_stator_v00.mat"); 
    N = length(data.Te_f_12);
    t = (0:N-1) / p.fs;
    idx = (t >= 2);
    t = t(idx);
    
    Te = data.Te_f_12(idx, 1);
    ia = data.isa_f_12(idx, 1);
    ib = data.isb_f_12(idx, 1);
    ic = data.isc_f_12(idx, 1);

    % Complex Space Vector (CSV)
    is = sqrt(2/3) * (ia + ib * exp(1i*2*pi/3) + ic * exp(1i*4*pi/3));
    
    % Paper-compliant Symmetrical Components (Hilbert + Moving RMS)
    [i_zero, i_pos, i_neg] = compute_hilbert_symmetrical_components(ia, ib, ic, p.fs, p.f0);
    
    signal = struct('t', t, 'is', is, 'i_zero', i_zero, 'i_pos', i_pos, 'i_neg', i_neg, 'Te', Te);
end

function signal = load_healthy_sample(p)
    data = load("healthy_v00.mat"); 
    N = length(data.Te_h);
    t = (0:N-1) / p.fs;
    idx = (t >= 2);
    t = t(idx);
    
    Te = data.Te_h(idx, 1);
    ia = data.isa_h(idx, 1);
    ib = data.isb_h(idx, 1);
    ic = data.isc_h(idx, 1);

    % Complex Space Vector (CSV)
    is = sqrt(2/3) * (ia + ib * exp(1i*2*pi/3) + ic * exp(1i*4*pi/3));
    
    % Paper-compliant Symmetrical Components (Hilbert + Moving RMS)
    [i_zero, i_pos, i_neg] = compute_hilbert_symmetrical_components(ia, ib, ic, p.fs, p.f0);
    
    signal = struct('t', t, 'is', is, 'i_zero', i_zero, 'i_pos', i_pos, 'i_neg', i_neg, 'Te', Te);
end

%% HILBERT + MOVING RMS SYMMETRICAL COMPONENTS
function [i_zero, i_pos, i_neg] = compute_hilbert_symmetrical_components(ia, ib, ic, fs, f0)
    theta_a = angle(hilbert(ia));
    theta_b = angle(hilbert(ib));
    theta_c = angle(hilbert(ic));

    W = round(fs / f0);
    Ia_rms = sqrt(movmean(ia.^2, W));
    Ib_rms = sqrt(movmean(ib.^2, W));
    Ic_rms = sqrt(movmean(ic.^2, W));

    I_a_vec = Ia_rms .* exp(1i * theta_a);
    I_b_vec = Ib_rms .* exp(1i * theta_b);
    I_c_vec = Ic_rms .* exp(1i * theta_c);

    a = exp(1i * 2 * pi / 3);
    i_zero = (1/3) * (I_a_vec + I_b_vec + I_c_vec);
    i_pos  = (1/3) * (I_a_vec + a * I_b_vec + a^2 * I_c_vec);
    i_neg  = (1/3) * (I_a_vec + a^2 * I_b_vec + a * I_c_vec);
end

%% TIME DOMAIN FEATURES
function val = compute_peak(x), val = max(abs(x)); end
function val = compute_rms(x), val = sqrt(mean(abs(x).^2)); end
function val = compute_kurtosis(x), e = abs(x(:)); val = mean(((e - mean(e)) / std(e, 1)).^4); end
function val = compute_crest_factor(x), val = compute_peak(x) / compute_rms(x); end
function val = compute_clearance_factor(x), e = abs(x(:)); val = max(e) / (mean(sqrt(e))^2); end
function val = compute_impulse_factor(x), e = abs(x(:)); val = max(e) / mean(e); end
function val = compute_shape_factor(x), e = abs(x(:)); val = compute_rms(x) / mean(e); end
function val = compute_skewness(x), e = abs(x(:)); val = mean(((e - mean(e)) / std(e, 1)).^3); end
function val = compute_square_mean_root(x), e = abs(x(:)); val = (mean(sqrt(e)))^2; end
function val = compute_5th_norm_moment(x), e = abs(x(:)); val = mean(((e - mean(e)) / std(e, 1)).^5); end
function val = compute_6th_norm_moment(x), e = abs(x(:)); val = mean(((e - mean(e)) / std(e, 1)).^6); end
function val = compute_mean(x), val = mean(abs(x)); end
function val = compute_shape_factor_2(x), val = compute_square_mean_root(x) / (compute_mean(x) + eps); end
function val = compute_peak_to_peak(x), e = abs(x(:)); val = max(e) - min(e); end
function val = compute_kurtosis_factor(x), val = compute_kurtosis(x) / (compute_rms(x)^4); end
function val = compute_std(x), val = std(abs(x(:)), 1); end
function val = compute_smoothness(x), val = 1 - (1 / (1 + var(abs(x(:)), 1))); end
function val = compute_uniformity(x), val = 1 - (compute_std(x) / (compute_mean(x) + eps)); end
function val = compute_neg_log_likelihood(x)
    e = abs(x(:)); mu = mean(e); sigma = std(e, 1);
    val = mean((1 / (sigma * sqrt(2*pi))) * exp(-((e - mu).^2) / (2 * sigma^2)));
end
function val = compute_ac_dc_ratio(x), e = abs(x(:)); val = std(e, 1) / abs(mean(e)); end
function val = compute_tkeo_mean(x)
    x = x(:); psi = abs(x(2:end-1)).^2 - real(x(1:end-2) .* conj(x(3:end)));
    val = mean(psi);
end
function val = compute_tkeo_kurtosis(x)
    x = x(:); psi = abs(x(2:end-1)).^2 - real(x(1:end-2) .* conj(x(3:end)));
    val = mean(((psi - mean(psi)) / std(psi, 1)).^4);
end

%% FREQUENCY DOMAIN FEATURES
function [X_mag, f] = get_complex_fft_spectrum(x, fs)
    x = x(:); N = length(x);
    X_shift = fftshift(fft(x)) / N;
    X_mag = abs(X_shift);
    f = (-floor(N/2) : ceil(N/2)-1)' * (fs / N);
end
function val = compute_spectral_centroid(x, fs)
    [X_mag, f] = get_complex_fft_spectrum(x, fs); val = sum(f .* X_mag) / sum(X_mag);
end
function val = compute_spectral_spread(x, fs)
    [X_mag, f] = get_complex_fft_spectrum(x, fs); f_sc = compute_spectral_centroid(x, fs);
    val = sqrt(sum(((f - f_sc).^2) .* X_mag) / sum(X_mag));
end
function val = compute_spectral_kurtosis(x, fs)
    [X_mag, ~] = get_complex_fft_spectrum(x, fs);
    val = mean(((X_mag - mean(X_mag)) / std(X_mag, 1)).^4);
end
function val = compute_target_band_energy(x, fs, f_lower, f_upper)
    [X_mag, f] = get_complex_fft_spectrum(x, fs);
    band_idx = (abs(f) >= f_lower) & (abs(f) <= f_upper);
    val = sum(X_mag(band_idx).^2) / sum(X_mag.^2);
end
function val = compute_harmonic_distortion(x, fs, f0, H)
    [X_mag, f] = get_complex_fft_spectrum(x, fs);
    [~, idx_pos_f0] = min(abs(f - f0)); mag_pos_f0 = X_mag(idx_pos_f0);
    mag_harmonics_sq = 0;
    for h = 2:H
        [~, idx_pos_h] = min(abs(f - h * f0));
        [~, idx_neg_h] = min(abs(f + h * f0));
        mag_harmonics_sq = mag_harmonics_sq + X_mag(idx_pos_h)^2 + X_mag(idx_neg_h)^2;
    end
    val = sqrt(mag_harmonics_sq) / mag_pos_f0;
end
function val = compute_spectral_flatness(x, fs)
    [X_mag, ~] = get_complex_fft_spectrum(x, fs);
    val = exp(mean(log(X_mag + 1e-12))) / mean(X_mag);
end
function val = compute_kurtogram_peak(x, fs)
    K = kurtogram(abs(x(:)), fs); val = max(K(:));
end

%% CROSS-COMPONENT & GEOMETRIC FEATURES
function val = compute_cuf(i_pos, i_neg)
    val = (sqrt(mean(abs(i_neg(:)).^2)) / sqrt(mean(abs(i_pos(:)).^2))) * 100;
end
function [a, b] = fit_ellipse_axes(i_qds)
    Z = [real(i_qds(:)), imag(i_qds(:))];
    [~, S, ~] = svd(Z - mean(Z, 1), 0);
    sigmas = diag(S) * sqrt(2 / length(i_qds));
    a = max(sigmas); b = min(sigmas);
end
function val = compute_eccentricity(i_qds)
    [a, b] = fit_ellipse_axes(i_qds); val = sqrt(1 - (b / a)^2);
end
function val = compute_orbit_thickness(i_qds)
    r = abs(i_qds(:)); val = std(r, 1) / mean(r);
end
function val = compute_enclosed_area(i_qds)
    [a, b] = fit_ellipse_axes(i_qds); val = pi * a * b;
end

%% NONLINEAR DYNAMICS & COMPLEXITY FEATURES
function val = compute_sample_entropy(x, m, r)
    e = abs(x(:)); N = length(e);
    if nargin < 2 || isempty(m), m = 2; end
    if nargin < 3 || isempty(r), r = 0.2 * std(e, 1); end
    Xm  = zeros(N - m, m); Xm1 = zeros(N - m, m + 1);
    for k = 1:m
        Xm(:, k)  = e(k : N - m + k - 1);
        Xm1(:, k) = e(k : N - m + k - 1);
    end
    Xm1(:, m+1) = e(m + 1 : N);
    B = 0; A = 0; N_m = N - m;
    for i = 1 : N_m - 1
        dist_m  = max(abs(Xm(i+1:end, :)  - Xm(i, :)),  [], 2);
        dist_m1 = max(abs(Xm1(i+1:end, :) - Xm1(i, :)), [], 2);
        B = B + sum(dist_m < r); A = A + sum(dist_m1 < r);
    end
    if A == 0 || B == 0, val = 0; else, val = -log(A / B); end
end

function val = compute_higuchi_fractal_dim(x, kmax)
    e = abs(x(:)); N = length(e);
    if nargin < 2 || isempty(kmax), kmax = 10; end
    L = zeros(kmax, 1);
    for k = 1:kmax
        L_m = zeros(k, 1);
        for m = 1:k
            idx = m:k:N; N_m = length(idx);
            if N_m > 1
                L_m(m) = (sum(abs(diff(e(idx)))) * (N - 1) / (floor((N - m) / k) * k)) / k;
            end
        end
        L(k) = mean(L_m);
    end
    p = polyfit(log(1 ./ (1:kmax)'), log(L), 1); val = p(1);
end

function val = compute_lempel_ziv_complexity(x)
    e = abs(x(:)); N = length(e); s = e >= median(e);
    c = 1; l = 1; i = 1; k = 1; k_max = 1;
    while (i + k) <= N
        if s(i + k) == s(l + k)
            k = k + 1;
            if (l + k) > i
                i = i + k_max; c = c + 1; l = 1; k = 1; k_max = 1;
            end
        else
            if k > k_max, k_max = k; end
            l = l + 1;
            if l == i
                i = i + k_max; c = c + 1; l = 1; k = 1; k_max = 1;
            else, k = 1; end
        end
    end
    val = c / (N / log2(N));
end

function val = compute_hilbert_envelope_entropy(x)
    e = abs(x(:)); E_total = sum(e.^2);
    if E_total == 0, val = 0; return; end
    p = (e.^2) / E_total; p = p(p > 0);
    val = -sum(p .* log2(p));
end

%% RESULT REORGANIZATION
function feature_struct_array = build_feature_structs(feature_names, sample_feature_evolution, window_times_all, db)
    n_features = length(feature_names);
    n_samples  = length(db);

    feature_struct_array = repmat(struct('name', '', 'samples', []), 1, n_features);

    for j = 1:n_features
        sample_data = repmat(struct('label', '', 'time', [], 'values', []), 1, n_samples);
        for k = 1:n_samples
            sample_data(k).label  = db(k).label;
            sample_data(k).time   = window_times_all{k};
            sample_data(k).values = sample_feature_evolution{k}(:, j);
        end
        feature_struct_array(j).name    = feature_names{j};
        feature_struct_array(j).samples = sample_data;
    end
end

%% RELATIVE VARIATION
function relvar_table = compute_relative_variation(feature_struct_array)
    n_features = length(feature_struct_array);

    names        = cell(n_features, 1);
    var_healthy  = zeros(n_features, 1);
    var_faulty   = zeros(n_features, 1);
    rel_variance = zeros(n_features, 1);

    for j = 1:n_features
        fj = feature_struct_array(j);
        labels = {fj.samples.label};

        idx_healthy = find(strcmpi(labels, 'Healthy'), 1);
        idx_faulty  = find(strcmpi(labels, 'Faulty'), 1);

        v_h = var(fj.samples(idx_healthy).values);
        v_f = var(fj.samples(idx_faulty).values);

        names{j}        = fj.name;
        var_healthy(j)  = v_h;
        var_faulty(j)   = v_f;
        rel_variance(j) = (v_f - v_h) / (v_h + eps) * 100;
    end

    relvar_table = table(names, var_healthy, var_faulty, rel_variance, ...
        'VariableNames', {'Feature', 'Var_Healthy', 'Var_Faulty', 'RelativeVariance_pct'});

    % Sort descending by absolute magnitude of variance change
    [~, sort_idx] = sort(abs(rel_variance), 'descend');
    relvar_table  = relvar_table(sort_idx, :);
end

%% SINGLE-FEATURE FIGURE BUILDER
function fig = plot_single_feature(fj, visibility)
    if nargin < 2 || isempty(visibility), visibility = 'on'; end

    fig = figure('Name', fj.name, 'NumberTitle', 'off', 'Visible', visibility);
    hold on; grid on; box on;

    legend_labels = cell(1, length(fj.samples));
    for k = 1:length(fj.samples)
        plot(fj.samples(k).time, fj.samples(k).values, '-o', ...
             'LineWidth', 1.2, 'MarkerSize', 3);
        legend_labels{k} = fj.samples(k).label;
    end

    title(sprintf('Evolution of %s', fj.name), 'Interpreter', 'none');
    xlabel('Time (s)');
    ylabel(fj.name, 'Interpreter', 'none');
    legend(legend_labels, 'Location', 'best');
    hold off;
end

%% PLOTTING: ON SCREEN
function plot_feature_evolutions(feature_struct_array)
    n_features = length(feature_struct_array);
    for j = 1:n_features
        plot_single_feature(feature_struct_array(j), 'on');
    end
end

%% EXPORT PDF REPORT
function export_figures_to_pdf(feature_struct_array, output_filepath)
    if nargin < 2 || isempty(output_filepath)
        script_dir = fileparts(mfilename('fullpath'));
        if isempty(script_dir), script_dir = pwd; end
        output_filepath = fullfile(script_dir, 'feature_evolution_report.pdf');
    end

    if exist(output_filepath, 'file')
        delete(output_filepath);
    end

    n_features = length(feature_struct_array);
    [out_dir, out_name] = fileparts(output_filepath);
    pages_dir = fullfile(out_dir, [out_name '_pages']);

    append_supported = true;

    for j = 1:n_features
        fj  = feature_struct_array(j);
        fig = plot_single_feature(fj, 'off');

        if append_supported
            try
                exportgraphics(fig, output_filepath, 'Append', true, 'ContentType', 'vector');
            catch
                append_supported = false;
                if exist(output_filepath, 'file'), delete(output_filepath); end
                if ~exist(pages_dir, 'dir'), mkdir(pages_dir); end
                page_name = sprintf('%02d_%s.pdf', j, matlab.lang.makeValidName(fj.name));
                exportgraphics(fig, fullfile(pages_dir, page_name), 'ContentType', 'vector');
            end
        else
            page_name = sprintf('%02d_%s.pdf', j, matlab.lang.makeValidName(fj.name));
            exportgraphics(fig, fullfile(pages_dir, page_name), 'ContentType', 'vector');
        end

        close(fig);
    end

    if append_supported
        fprintf('PDF report exported: %d figures saved to:\n%s\n', n_features, output_filepath);
    else
        fprintf(['This MATLAB release does not support multi-page PDF append.\n' ...
                 '%d single-page PDFs were saved instead in:\n%s\n'], ...
                 n_features, pages_dir);
    end
end

%% PLOTTING: TOP N FEATURES BAR CHART (For Mentor Presentation Slides)
function plot_top_relative_variance_barchart(relvar_table, top_N)
    if nargin < 2 || isempty(top_N), top_N = 25; end
    
    % Select top N features by absolute magnitude
    T_sub = relvar_table(1:min(top_N, height(relvar_table)), :);
    
    % Flip rows so highest feature appears at the TOP of the horizontal bar plot
    T_plot = flipud(T_sub);
    n_items = height(T_plot);

    figure('Name', sprintf('Top %d Relative Variance Features', top_N), 'NumberTitle', 'off', ...
           'Units', 'pixels', 'Position', [150, 150, 950, 700]);

    % Plot bars directly with explicit numeric indices to prevent categorical label inversion
    hBar = barh(1:n_items, T_plot.RelativeVariance_pct);
    hBar.FaceColor = [0.18 0.45 0.71];

    yticks(1:n_items);
    yticklabels(T_plot.Feature);

    xlabel('Relative Variance Change (%)', 'FontSize', 12, 'FontWeight', 'bold');
    ylabel('Feature Name', 'FontSize', 12, 'FontWeight', 'bold');
    title(sprintf('Top %d Diagnostic Features Sorted by Absolute Sensitivity (Faulty vs. Healthy)', top_N), ...
          'FontSize', 13, 'FontWeight', 'bold', 'Interpreter', 'none');

    grid on; box on;

    % Vertical reference line at 0%
    xline(0, 'k--', 'LineWidth', 1.2);

    ax = gca;
    ax.YAxis.FontSize = 9.5;
    ax.YAxis.TickLabelInterpreter = 'none';
    ylim([0.5, n_items + 0.5]);
    xlim padded;
end

%% PLOTTING: ALL 117 FEATURES BAR CHART (Full Reference Plot)
function plot_relative_variance_barchart(relvar_table)
    % Flip rows so highest feature appears at the TOP of the plot
    T_plot = flipud(relvar_table);
    n_features = height(T_plot);

    % Tall figure window to accommodate all 117 tick labels cleanly
    figure('Name', 'Relative Variance (All Features)', 'NumberTitle', 'off', ...
           'Units', 'pixels', 'Position', [100, 50, 1200, 2400]);

    % Plot bars directly with numeric indices to preserve exact row matching
    hBar = barh(1:n_features, T_plot.RelativeVariance_pct);
    hBar.FaceColor = [0.18 0.45 0.71];

    yticks(1:n_features);
    yticklabels(T_plot.Feature);

    xlabel('Relative Variance Change (%)', 'FontSize', 11, 'FontWeight', 'bold');
    ylabel('Feature Name', 'FontSize', 11, 'FontWeight', 'bold');
    title(sprintf('Relative Variance Change Across All %d Features (Faulty vs. Healthy)', n_features), ...
          'FontSize', 12, 'FontWeight', 'bold', 'Interpreter', 'none');

    grid on; box on;

    % Vertical reference line at 0% to clearly display negative changes
    xline(0, 'k--', 'LineWidth', 1.2);

    ax = gca;
    ax.YAxis.FontSize = 6.5;
    ax.YAxis.TickLabelInterpreter = 'none';
    ylim([0.5, n_features + 0.5]);
    xlim padded;
end