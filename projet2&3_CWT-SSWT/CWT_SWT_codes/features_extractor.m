clc; clear; close all;

data = load('healthy_v00.mat'); 

% keep only the last 8s 
fs=5000;
N = 50001;
t = (0:N-1) / fs;
idx = (t >= 2);   
t = t(idx);
    
for i = 1:length(data)
    
    % Extraction des champs requis
    Te = data.Te_h(:, i);
    ia = data.isa_h(:, i);
    ib = data.isb_h(:, i);
    ic = data.isc_h(:, i);
    
    % keep only useful 8s
    Te = Te(idx);
    ia = ia(idx);
    ib = ib(idx);
    ic = ic(idx);
     

    A_pos = 5;      f_pos = +100;    % amplitude & frequency of +200 Hz component
    A_neg = 1;      f_neg = -100;    % amplitude & frequency of -200 Hz component

    % Complex space vector on the full timeline, then trim
    is = A_pos * exp(1i*2*pi*f_pos*t) + ...
              A_neg * exp(1i*2*pi*f_neg*t);


    
    % Compute the symmetrical components (Zero, Positive, and Negative sequences)
    a = exp(1i*2*pi/3);
    i_zero = (1/3) * (ia + ib + ic);
    i_pos  = (1/3) * (ia + a*ib + a^2*ic);
    i_neg  = (1/3) * (ia + a^2*ib + a*ic);
    
    %% compute and extract features:
    
    % static features related to ia, ib, ic:
    ia_features = extract_features(t, ia);
    ib_features = extract_features(t, ib);
    ic_features = extract_features(t, ic);
    
    % features related to sequences:
    seq_features = extract_sequence_features(i_zero, i_pos, i_neg);
    
    % static featuresn related to is:
    is_features = extract_space_vector_features(is);
    
    % obtain CWT results to extract features:
    results = Wavelet_Spectrum_Shifter(is, t);
    % extract feats from cwt analysis:
    cwt_features = extract_cwt_features(results);
    
    % extract features from the torque signal Te:
    torque_features = extract_torque_features(Te, fs, 50.0); 
    
    %% --- 4. Unpack Every Single Elementary Feature Individually ---
    
    % Phase A Features (Columns 1-7)
    f_ia_RMS          = ia_features.RMS;
    f_ia_Peak         = ia_features.Peak;
    f_ia_CrestFactor  = ia_features.CrestFactor;
    f_ia_Kurtosis     = ia_features.Kurtosis;
    f_ia_Skewness     = ia_features.Skewness;
    f_ia_ShapeFactor  = ia_features.ShapeFactor;
    f_ia_MarginFactor = ia_features.MarginFactor;
    
    % Phase B Features (Columns 8-14)
    f_ib_RMS          = ib_features.RMS;
    f_ib_Peak         = ib_features.Peak;
    f_ib_CrestFactor  = ib_features.CrestFactor;
    f_ib_Kurtosis     = ib_features.Kurtosis;
    f_ib_Skewness     = ib_features.Skewness;
    f_ib_ShapeFactor  = ib_features.ShapeFactor;
    f_ib_MarginFactor = ib_features.MarginFactor;
    
    % Phase C Features (Columns 15-21)
    f_ic_RMS          = ic_features.RMS;
    f_ic_Peak         = ic_features.Peak;
    f_ic_CrestFactor  = ic_features.CrestFactor;
    f_ic_Kurtosis     = ic_features.Kurtosis;
    f_ic_Skewness     = ic_features.Skewness;
    f_ic_ShapeFactor  = ic_features.ShapeFactor;
    f_ic_MarginFactor = ic_features.MarginFactor;
    
    % Symmetrical Sequence Features (Columns 22-37)
    f_seq_neg_std                   = seq_features.neg_std;
    f_seq_zero_std                  = seq_features.zero_std;
    f_seq_Negative_Sequence_RMS     = seq_features.Negative_Sequence_RMS;
    f_seq_Zero_Sequence_RMS         = seq_features.Zero_Sequence_RMS;
    f_seq_Current_Unbalance_Factor  = seq_features.Current_Unbalance_Factor;
    f_seq_Zero_To_Positive_Ratio    = seq_features.Zero_To_Positive_Ratio;
    f_seq_Negative_Sequence_Kurtosis= seq_features.Negative_Sequence_Kurtosis;
    f_seq_Zero_Sequence_Kurtosis    = seq_features.Zero_Sequence_Kurtosis;
    f_seq_Phase_Jitter_Std          = seq_features.Phase_Jitter_Std;
    f_seq_Phase_Diff_Neg_Pos        = seq_features.Phase_Diff_Neg_Pos;
    f_seq_Phase_Diff_Zero_Neg       = seq_features.Phase_Diff_Zero_Neg;
    f_seq_Negative_Crest_Factor     = seq_features.Negative_Crest_Factor;
    f_seq_Zero_Crest_Factor         = seq_features.Zero_Crest_Factor;
    f_seq_Positive_Sequence_Variance= seq_features.Positive_Sequence_Variance;
    f_seq_Negative_Sequence_Variance= seq_features.Negative_Sequence_Variance;
    f_seq_Zero_Sequence_Variance    = seq_features.Zero_Sequence_Variance;
    % Note: If using Signal Processing Toolbox for THD, add those here as well.
    
    % Complex Space Vector Features (Columns 38-45)
    f_csv_Mean           = is_features.Mean;
    f_csv_RippleVariance = is_features.RippleVariance;
    f_csv_P2PMargin      = is_features.P2PMargin;
    f_csv_Eccentricity   = is_features.Eccentricity;
    f_csv_Kurtosis       = is_features.Kurtosis;
    f_csv_Skewness       = is_features.Skewness;
    f_csv_ShapeFactor    = is_features.ShapeFactor;
    f_csv_RMS            = is_features.RMS;
    
    % Wavelet (CWT) Scalogram Features (Columns 46-53)
    f_cwt_total_energy             = cwt_features.total_energy;
    f_cwt_positive_energy          = cwt_features.positive_energy;
    f_cwt_negative_energy          = cwt_features.negative_energy;
    f_cwt_energy_ratio             = cwt_features.energy_ratio;
    f_cwt_positive_entropy         = cwt_features.positive_entropy;
    f_cwt_negative_entropy         = cwt_features.negative_entropy;
    f_cwt_spectral_centroid        = cwt_features.spectral_centroid;
    f_cwt_target_band_energy_ratio = cwt_features.target_band_energy_ratio;
    
    % Mechanical Torque Features (Columns 54-63)
    f_torque_Mean         = torque_features.Torque_Mean;
    f_torque_std          = torque_features.Torque_std;
    f_torque_RMS_AC       = torque_features.Torque_RMS_AC;
    f_torque_Peak2Peak    = torque_features.Torque_Peak2Peak;
    f_torque_Variance     = torque_features.Torque_Variance;
    f_torque_Kurtosis     = torque_features.Torque_Kurtosis;
    f_torque_AC_DC_Ratio  = torque_features.Torque_AC_DC_Ratio;
    f_torque_Crest_Factor = torque_features.Torque_Crest_Factor;
    f_torque_Ripple_2xGrid= torque_features.Torque_Ripple_2xGrid;
    f_torque_THD          = torque_features.Torque_THD;
    
    %%  Assemble Flat Raw Row into Matrix ---
    Feature_Matrix(i, :) = [ ...
        f_ia_RMS, f_ia_Peak, f_ia_CrestFactor, f_ia_Kurtosis, f_ia_Skewness, f_ia_ShapeFactor, f_ia_MarginFactor, ...
        f_ib_RMS, f_ib_Peak, f_ib_CrestFactor, f_ib_Kurtosis, f_ib_Skewness, f_ib_ShapeFactor, f_ib_MarginFactor, ...
        f_ic_RMS, f_ic_Peak, f_ic_CrestFactor, f_ic_Kurtosis, f_ic_Skewness, f_ic_ShapeFactor, f_ic_MarginFactor, ...
        f_seq_neg_std, f_seq_zero_std, f_seq_Negative_Sequence_RMS, f_seq_Zero_Sequence_RMS, ...
        f_seq_Current_Unbalance_Factor, f_seq_Zero_To_Positive_Ratio, f_seq_Negative_Sequence_Kurtosis, ...
        f_seq_Zero_Sequence_Kurtosis, f_seq_Phase_Jitter_Std, f_seq_Phase_Diff_Neg_Pos, f_seq_Phase_Diff_Zero_Neg, ...
        f_seq_Negative_Crest_Factor, f_seq_Zero_Crest_Factor, f_seq_Positive_Sequence_Variance, ...
        f_seq_Negative_Sequence_Variance, f_seq_Zero_Sequence_Variance, ...
        f_csv_Mean, f_csv_RippleVariance, f_csv_P2PMargin, f_csv_Eccentricity, f_csv_Kurtosis, f_csv_Skewness, f_csv_ShapeFactor, f_csv_RMS, ...
        f_cwt_total_energy, f_cwt_positive_energy, f_cwt_negative_energy, f_cwt_energy_ratio, f_cwt_positive_entropy, f_cwt_negative_entropy, f_cwt_spectral_centroid, f_cwt_target_band_energy_ratio, ...
        f_torque_Mean, f_torque_std, f_torque_RMS_AC, f_torque_Peak2Peak, f_torque_Variance, f_torque_Kurtosis, f_torque_AC_DC_Ratio, f_torque_Crest_Factor, f_torque_Ripple_2xGrid, f_torque_THD ...
    ];

%     if mod(i, 10) == 0 || i == num_samples
%         fprintf('Successfully processed matrix row %d/%d\n', i, num_samples);
%     end
end

% Save Complete Array Data ---

% Get the folder where this specific script/function is saved
[current_folder, ~, ~] = fileparts(mfilename('fullpath'));
    
% Combine the folder path with your desired filename
forced_filepath = fullfile(current_folder, 'extracted_elementary_features.mat');
    
% Save the results struct
save(forced_filepath, 'Feature_Matrix');
    
fprintf('\nDone! Feature_Matrix saved successfully as an isolated 2D numeric grid.\n');
fprintf('\n[SAVED] results written to %s\n', forced_filepath);


%% features related to Ia Ib Ic :
function features = extract_features(t, signal)

    signal = signal(:); 

    % --- 1. Call the Feature Functions ---
    features.RMS          = compute_RMS(signal);
    features.Peak         = find_peak(signal);
    features.CrestFactor  = compute_Crest_Factor(signal);
    features.Kurtosis     = compute_Kurtosis(signal);
    features.Skewness     = compute_Skewness(signal);
    features.ShapeFactor  = compute_shape_Factor(signal);
    features.MarginFactor = compute_Margin_Factor(signal);
end

% =========================================================================
% INDIVIDUAL FEATURE FUNCTIONS (Local Functions)
% =========================================================================

function rms_val = compute_RMS(x)
    % Root Mean Square: Square root of the arithmetic mean of the squares
    rms_val = sqrt(mean(x.^2));
end

function peak_val = find_peak(x)
    % Maximum absolute peak value of the signal
    peak_val = max(abs(x));
end

function cf = compute_Crest_Factor(x)
    % Crest Factor: Peak amplitude divided by the RMS value
    % High values indicate extreme impulses/spikes in the signal
    peak_val = find_peak(x);
    rms_val = compute_RMS(x);
    cf = peak_val / rms_val;
end

function kurt_val = compute_Kurtosis(x)
    % Kurtosis: Measures the "tailedness" of the probability distribution
    % Uses MATLAB's built-in function. (3 = perfectly normal/Gaussian distribution)
    kurt_val = kurtosis(x);
end

function skew_val = compute_Skewness(x)
    % Skewness: Measures the asymmetry of the data distribution around its mean
    skew_val = skewness(x);
end

function sf = compute_shape_Factor(x)
    % Shape Factor: RMS value divided by the Mean Absolute value
    % Quantifies the structural shape of the wave profile
    rms_val = compute_RMS(x);
    mean_abs = mean(abs(x));
    sf = rms_val / mean_abs;
end

function mf = compute_Margin_Factor(x)
    % Margin Factor: Peak amplitude divided by the Square of the Mean of the Square Roots
    % Highly sensitive to changes in early-stage component degradation
    peak_val = find_peak(x);
    root_mean_square_root = (mean(sqrt(abs(x)))).^2;
    mf = peak_val / root_mean_square_root;
end

%% featuresrelated to zero/ pos/ neg Seq
function seq_features = extract_sequence_features(i_zero, i_pos, i_neg)

    % --- Your Original Features ---
    seq_features.neg_std                   = std(abs(i_neg));
    seq_features.zero_std                  = std(abs(i_zero));
    seq_features.Negative_Sequence_RMS      = compute_negative_Sequence_RMS(i_neg);
    seq_features.Zero_Sequence_RMS          = compute_zero_Sequence_RMS(i_zero);
    seq_features.Current_Unbalance_Factor   = compute_Current_Unbalance_Factor(i_pos, i_neg);
    seq_features.Zero_To_Positive_Ratio     = compute_Zero_To_Positive_Ratio(i_pos, i_zero);
    seq_features.Negative_Sequence_Kurtosis = compute_Negative_Sequence_Kurtosis(i_neg);
    seq_features.Zero_Sequence_Kurtosis     = compute_Zero_Sequence_Kurtosis(i_zero);
    
    % Inside extract_sequence_features:
    seq_features.inst_phase_diff = angle(i_neg) - angle(i_pos);
    seq_features.Phase_Jitter_Std = std(unwrap(seq_features.inst_phase_diff)) * (180/pi);

    % --- NEW: Phase Angular Features ---
    [phi_21, phi_02] = compute_Sequence_Phase_Diff(i_pos, i_neg, i_zero);
    seq_features.Phase_Diff_Neg_Pos         = phi_21;
    seq_features.Phase_Diff_Zero_Neg        = phi_02;

    % --- NEW: Crest Factor Features (Impulsive Energy) ---
    seq_features.Negative_Crest_Factor      = compute_Seq_Crest_Factor(i_neg);
    seq_features.Zero_Crest_Factor          = compute_Seq_Crest_Factor(i_zero);

    % --- NEW: Envelope Variance Features ---
    seq_features.Positive_Sequence_Variance = var(abs(i_pos));
    seq_features.Negative_Sequence_Variance = var(abs(i_neg));
    seq_features.Zero_Sequence_Variance     = var(abs(i_zero));

    % --- NEW: Harmonic Distortion Features (Requires Signal Processing Toolbox) ---
    seq_features.Negative_Sequence_THD      = thd(abs(i_neg));
    seq_features.Zero_Sequence_THD          = thd(abs(i_zero));
    
end

% =========================================================================
% INDIVIDUAL SEQUENCE FEATURE FUNCTIONS (Local Functions)
% =========================================================================

function rms_neg = compute_negative_Sequence_RMS(i_neg)
    % Root Mean Square of the complex negative sequence signal
    rms_neg = sqrt(mean(abs(i_neg).^2));
end

function rms_zero = compute_zero_Sequence_RMS(i_zero)
    % Root Mean Square of the complex zero sequence signal
    rms_zero = sqrt(mean(abs(i_zero).^2));
end

function cuf = compute_Current_Unbalance_Factor(i_pos, i_neg)
    % CUF (%) = (RMS of Negative Sequence / RMS of Positive Sequence) * 100
    % Standard industrial indicator for phase unbalance (e.g., stator faults)
    rms_pos = sqrt(mean(abs(i_pos).^2));
    rms_neg = sqrt(mean(abs(i_neg).^2));
    
    if rms_pos < 1e-6
        cuf = 0;
    else
        cuf = (rms_neg / rms_pos) * 100; 
    end
end

function zpr = compute_Zero_To_Positive_Ratio(i_pos, i_zero)
    % ZPR (%) = (RMS of Zero Sequence / RMS of Positive Sequence) * 100
    % Highly sensitive indicator for ground faults or neutral line leaks
    rms_pos = sqrt(mean(abs(i_pos).^2));
    rms_zero = sqrt(mean(abs(i_zero).^2));
    
    if rms_pos < 1e-6
        zpr = 0;
    else
        zpr = (rms_zero / rms_pos) * 100;
    end
end

function kurt_neg = compute_Negative_Sequence_Kurtosis(i_neg)
    % Kurtosis of the Negative Sequence magnitude profile
    % Useful to see if the unbalance is steady or impulsive/spiky
    kurt_neg = kurtosis(abs(i_neg));
end

function kurt_zero = compute_Zero_Sequence_Kurtosis(i_zero)
    % Kurtosis of the Zero Sequence magnitude profile
    % Spikes in ground leakage will yield high Kurtosis values here
    kurt_zero = kurtosis(abs(i_zero));
end

function [phi_21, phi_02] = compute_Sequence_Phase_Diff(i_pos, i_neg, i_zero)
    % Extracts mean circular phase angle configurations in degrees (-180 to 180)
    angle_pos  = mean(angle(i_pos)) * (180/pi);
    angle_neg  = mean(angle(i_neg)) * (180/pi);
    angle_zero = mean(angle(i_zero)) * (180/pi);
    
    % Compute phase discrepancies while managing circular wrapping bounds
    phi_21 = mod(angle_neg - angle_pos + 180, 360) - 180;
    phi_02 = mod(angle_zero - angle_neg + 180, 360) - 180;
end

function cf = compute_Seq_Crest_Factor(i_seq)
    % Compares peak dynamic amplitude directly against baseline RMS energy
    rms_val = sqrt(mean(abs(i_seq).^2));
    if rms_val < 1e-6
        cf = 1; % Crest factor of a zeroed signal defaults to unity baseline
    else
        cf = max(abs(i_seq)) / rms_val;
    end
end

%% features related to CSV
function CSV_Features = extract_space_vector_features(is)
    % Extracting the metrics
    CSV_Features.Mean           = compute_mean_CSV(is);
    CSV_Features.RippleVariance = compute_envelop_ripple_variance_CSV(is);
    CSV_Features.P2PMargin      = compute_envelop_P2P_Margin_CSV(is);
    CSV_Features.Eccentricity   = compute_Geometric_Eccentricity_CSV(is);
    CSV_Features.Kurtosis       = compute_kurtosis_CSV(is);
    CSV_Features.Skewness       = compute_skewness_CSV(is);
    CSV_Features.ShapeFactor    = compute_Shape_Factor_CSV(is);
    CSV_Features.RMS            = compute_RMS_CSV(is);

end

% =========================================================================
% SPECIALIZED COMPLEX SPACE VECTOR (CSV) FEATURE FUNCTIONS
% =========================================================================

function mean_csv = compute_mean_CSV(is)
    % Computes the mean of the Complex Space Vector envelope (magnitude)
    % Input: 'is' is the complex vector (i_alpha + j*i_beta)
    mean_csv = mean(abs(is));
end

function ripple_var = compute_envelop_ripple_variance_CSV(is)
    % Computes the variance of the envelope ripple.
    % In a healthy motor, the envelope (magnitude) is a constant radius (variance near 0).
    % Under a fault, the envelope ripples; this measures the severity of that ripple.
    env = abs(is);
    ripple_var = var(env);
end

function p2p_margin = compute_envelop_P2P_Margin_CSV(is)
    % Computes the Peak-to-Peak Margin Factor of the space vector envelope.
    % Formula: Peak-to-Peak (max - min) divided by the square of the mean of square roots.
    env = abs(is);
    p2p = max(env) - min(env);
    
    % Avoid division by zero if the signal is perfectly flat/zero
    root_mean_square_root = (mean(sqrt(env))).^2;
    if root_mean_square_root == 0
        p2p_margin = 0;
    else
        p2p_margin = p2p / root_mean_square_root;
    end
end

function eccentricity = compute_Geometric_Eccentricity_CSV(is)
    % Computes the geometric eccentricity of the space vector trajectory.
    % Healthy trajectory = perfect circle (eccentricity = 0).
    % Faulty trajectory = ellipse (eccentricity approaches 1).
    % Formula: e = sqrt(1 - (b^2 / a^2)) where 'a' is the semi-major axis and 'b' is the semi-minor axis.
    
    i_alpha = real(is);
    i_beta  = imag(is);
    
    % Use standard deviations of the components as an estimate for the axes
    std_alpha = std(i_alpha);
    std_beta  = std(i_beta);
    
    % Define major (a) and minor (b) based on which component fluctuates more
    a = max(std_alpha, std_beta);
    b = min(std_alpha, std_beta);
    
    if a == 0
        eccentricity = 0;
    else
        eccentricity = sqrt(1 - (b^2 / a^2));
    end
end

function kurt_csv = compute_kurtosis_CSV(is)
    % Computes the Kurtosis of the Space Vector magnitude profile.
    % Detects if the space vector is experiencing sharp, transient shocks.
    kurt_csv = kurtosis(abs(is));
end

function skew_csv = compute_skewness_CSV(is)
    % Computes the Skewness (asymmetry) of the Space Vector magnitude profile.
    skew_csv = skewness(abs(is));
end

function sf_csv = compute_Shape_Factor_CSV(is)
    % Computes the Shape Factor of the Space Vector envelope.
    % Formula: RMS of envelope / Mean of Absolute envelope
    env = abs(is);
    rms_env = sqrt(mean(env.^2));
    mean_abs_env = mean(env);
    
    if mean_abs_env == 0
        sf_csv = 0;
    else
        sf_csv = rms_env / mean_abs_env;
    end
end

function rms_csv = compute_RMS_CSV(is)
    % Computes the total Root Mean Square energy of the Complex Space Vector.
    % For a complex vector, this evaluates both real and imaginary energy components.
    rms_csv = sqrt(mean(abs(is).^2));
end


%% features extracted from CWT:
function feats = extract_cwt_features(results)

    % Extract variables from the cwt_info structure
    entry = results.cwt_info;
    W = entry.W;
    f = entry.f;
    
    feats.total_energy      = compute_total_Scalogram_Energy(W);
    feats.positive_energy   = compute_total_Positive_Scalogram_Energy(W);
    feats.negative_energy   = compute_total_Negative_Scalogram_Energy(W);
    feats.energy_ratio      = compute_Scalogram_Energy_Ratio(W);
    feats.positive_entropy  = compute_Positive_Shannon_Entropy(W);
    feats.negative_entropy  = compute_Negative_Shannon_Entropy(W);
    feats.spectral_centroid = compute_Spectral_Center_of_Gravity(W, f);
   
    
    feats.target_band_energy_ratio = compute_Target_Band_Energy_Ratio(W, entry.band_pos.idx, entry.band_neg.idx);
    feats.target_band_entropy_pos  = compute_Target_Band_Entropy(W, entry.band_pos.idx, 1);
    feats.target_band_entropy_neg  = compute_Target_Band_Entropy(W, entry.band_neg.idx, 2);

    
    feats.rec_energy_ratio  = compute_Reconstructed_Energy_Ratio(entry.x_rec_pos, entry.x_rec_neg);
    feats.rec_crest_factor_pos = compute_cwt_Crest_Factor(entry.x_rec_pos);
    feats.rec_crest_factor_neg = compute_cwt_Crest_Factor(entry.x_rec_neg);

    
    feats.fft_amp_ratio     = entry.amp_neg / (entry.amp_pos + eps);
    feats.fft_freq_asymmetry = entry.freq_pos - abs(entry.freq_neg);
    feats.pipeline_beta     = entry.beta;
    feats.pipeline_vpo      = entry.vpo;

end


function E_total = compute_total_Scalogram_Energy(W)
    E_total = sum(abs(W(:)).^2);
end

function E_pos = compute_total_Positive_Scalogram_Energy(W)
    W_pos = W(:,:,1);
    E_pos = sum(abs(W_pos(:)).^2);
end

function E_neg = compute_total_Negative_Scalogram_Energy(W)
    W_neg = W(:,:,2);
    E_neg = sum(abs(W_neg(:)).^2);
end

function ratio = compute_Scalogram_Energy_Ratio(W)
    E_pos = compute_total_Positive_Scalogram_Energy(W);
    E_neg = compute_total_Negative_Scalogram_Energy(W);
    ratio = E_neg / (E_pos + eps);
end

function H_pos = compute_Positive_Shannon_Entropy(W)
    W_pos = W(:,:,1);
    E = abs(W_pos(:)).^2;
    p = E / (sum(E) + eps);
    p = p(p > 0);
    H_pos = -sum(p .* log2(p));
end

function H_neg = compute_Negative_Shannon_Entropy(W)
    W_neg = W(:,:,2);
    E = abs(W_neg(:)).^2;
    p = E / (sum(E) + eps);
    p = p(p > 0);
    H_neg = -sum(p .* log2(p));
end

function fc = compute_Spectral_Center_of_Gravity(W, f)
    f = f(:);
    E_pos_per_scale = sum(abs(W(:,:,1)).^2, 2);
    E_neg_per_scale = sum(abs(W(:,:,2)).^2, 2);

    weighted_sum = sum(f .* E_pos_per_scale) + sum(-f .* E_neg_per_scale);
    total_energy = sum(E_pos_per_scale) + sum(E_neg_per_scale);

    fc = weighted_sum / (total_energy + eps);
end

function target_ratio = compute_Target_Band_Energy_Ratio(W, idx_pos, idx_neg)
    % Evaluates energy imbalances strictly inside your optimization bands
    W_target_pos = W(idx_pos, :, 1);
    W_target_neg = W(idx_neg, :, 2);
    
    E_target_pos = sum(abs(W_target_pos(:)).^2);
    E_target_neg = sum(abs(W_target_neg(:)).^2);
    
    target_ratio = E_target_neg / (E_target_pos + eps);
end

function H_band = compute_Target_Band_Entropy(W, idx, page)
    % Tracks transient complexity exclusively inside the filtered frequency band
    W_target = W(idx, :, page);
    E = abs(W_target(:)).^2;
    p = E / (sum(E) + eps);
    p = p(p > 0);
    H_band = -sum(p .* log2(p));
end

function rec_ratio = compute_Reconstructed_Energy_Ratio(x_pos, x_neg)
    % Captures the energy breakdown of the forward vs backward time-domain waveforms
    E_pos = sum(abs(x_pos(:)).^2);
    E_neg = sum(abs(x_neg(:)).^2);
    rec_ratio = E_neg / (E_pos + eps);
end

function cf = compute_cwt_Crest_Factor(x_rec)
    % Measures the impulsiveness/spikiness of the localized component
    rms_val = sqrt(mean(abs(x_rec).^2));
    if rms_val < 1e-6
        cf = 1; 
    else
        cf = max(abs(x_rec)) / rms_val;
    end
end


%% features related to the torque signal:
function torque_feats = extract_torque_features(torque, fs, grid_freq)
    
    % Ensure consistent column vector orientation
    torque = torque(:);
    
    % Decompose signal into operational DC (load) and dynamic AC (ripples)
    t_dc = mean(torque);
    t_ac = torque - t_dc;
    
    
    torque_feats.Torque_Mean         = t_dc;
    torque_feats.Torque_std          = std(torque);
    torque_feats.Torque_RMS_AC       = compute_Torque_RMS_AC(t_ac);
    torque_feats.Torque_Peak2Peak    = compute_Torque_Peak2Peak(torque);
    torque_feats.Torque_Variance     = compute_Torque_Variance(torque);
    torque_feats.Torque_Kurtosis     = compute_Torque_Kurtosis(torque);
    torque_feats.Torque_AC_DC_Ratio  = compute_Torque_AC_DC_Ratio(t_dc, t_ac);
    torque_feats.Torque_Crest_Factor = compute_Torque_Crest_Factor(t_ac);
    torque_feats.Torque_Ripple_2xGrid = compute_Torque_2xGrid_Ripple(t_ac, fs, grid_freq);
    torque_feats.Torque_THD          = compute_Torque_THD(t_ac);

end


function rms_ac = compute_Torque_RMS_AC(t_ac)
    % Computes the Root Mean Square of the isolated alternating components
    rms_ac = sqrt(mean(t_ac.^2));
end

function p2p = compute_Torque_Peak2Peak(torque)
    % Measures the total absolute range/span of the twisting force profile
    p2p = max(torque) - min(torque);
end

function t_var = compute_Torque_Variance(torque)
    % Captures the dispersion/instability of the torque over the window
    t_var = var(torque);
end

function t_kurt = compute_Torque_Kurtosis(torque)
    % Measures the heavy-tailed, spiky nature of severe physical impacts
    t_kurt = kurtosis(torque);
end

function ac_dc_ratio = compute_Torque_AC_DC_Ratio(t_dc, t_ac)
    % Normalizes structural ripple variations against the current base load
    rms_ac = sqrt(mean(t_ac.^2));
    ac_dc_ratio = rms_ac / (abs(t_dc) + eps);
end

function cf = compute_Torque_Crest_Factor(t_ac)
    % Tracks extreme transient spikes relative to typical steady ripple levels
    rms_ac = sqrt(mean(t_ac.^2));
    if rms_ac < 1e-6
        cf = 1; 
    else
        cf = max(abs(t_ac)) / rms_ac;
    end
end

function ripple_2x = compute_Torque_2xGrid_Ripple(t_ac, fs, grid_freq)
    % Extracts the electrical unbalance magnitude signature (2 * Grid Frequency)
    N = length(t_ac);
    Y = fft(t_ac) / N;
    f_axis = (0:N-1) * (fs / N);
    
    % Keep only the single-sided positive spectrum
    amp_spectrum = 2 * abs(Y(1:floor(N/2)+1));
    f_axis = f_axis(1:floor(N/2)+1);
    
    % Track targeted 2x grid marker within a narrow +/- 2 Hz band
    target_f = 2 * grid_freq;
    idx_target = (f_axis >= target_f - 2) & (f_axis <= target_f + 2);
    
    if any(idx_target)
        ripple_2x = max(amp_spectrum(idx_target));
    else
        ripple_2x = 0;
    end
end

function t_thd = compute_Torque_THD(t_ac)
    % Quantifies global harmonic contamination due to arcs or structural degradation
    % (Requires Signal Processing Toolbox)
    t_thd = thd(t_ac);
end