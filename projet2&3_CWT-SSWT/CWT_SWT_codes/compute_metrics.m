
%% features related to Ia Ib Ic :
function features = extract_features(t, signal)
% MAIN FUNCTION: Extracts all 7 features from a given signal
% Inputs:
%   t      - Time vector
%   signal - Signal values (1D vector)
% Output:
%   features - A structure containing all calculated metrics

    % Ensure the signal is a column vector for consistent matrix math
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

%% Feature related to positive/zero/negative Sequences :

function seq_features = extract_sequence_features(t, ia, ib, ic)
% MAIN FUNCTION: Computes symmetrical components and extracts sequence features
% Inputs:
%   t          - Time vector
%   ia, ib, ic - Three-phase current vectors (must be same length)
% Output:
%   seq_features - A structure containing all calculated metrics

    % Ensure inputs are column vectors
    ia = ia(:); ib = ib(:); ic = ic(:);
    
    % --- 1. Symmetrical Components Transformation (Fortescue) ---
    % Complex operator a = exp(j * 120 degrees)
    a = exp(1j * 2 * pi / 3);
    
    % Compute instantaneous symmetrical components sequence vectors
    % i_zero = 1/3 * (ia + ib + ic)
    % i_pos  = 1/3 * (ia + a*ib + a^2*ic)
    % i_neg  = 1/3 * (ia + a^2*ib + a*ic)
    i_zero = (1/3) * (ia + ib + ic);
    i_pos  = (1/3) * (ia + a * ib + (a^2) * ic);
    i_neg  = (1/3) * (ia + (a^2) * ib + a * ic);

    % --- 2. Call the Sequence Feature Functions ---
    seq_features.Negative_Sequence_RMS      = compute_negative_Sequence_RMS(i_neg);
    seq_features.Zero_Sequence_RMS          = compute_zero_Sequence_RMS(i_zero);
    seq_features.Current_Unbalance_Factor   = compute_Current_Unbalance_Factor(i_pos, i_neg);
    seq_features.Zero_To_Positive_Ratio     = compute_Zero_To_Positive_Ratio(i_pos, i_zero);
    seq_features.Negative_Sequence_Kurtosis = compute_Negative_Sequence_Kurtosis(i_neg);
    seq_features.Zero_Sequence_Kurtosis     = compute_Zero_Sequence_Kurtosis(i_zero);

    % --- 3. Display Results in Command Window ---
    fprintf('\n====== EXTRACTED SEQUENCE FEATURES ======\n');
    disp(struct2table(seq_features));
    fprintf('=========================================\n');
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
    
    if rms_pos == 0
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
    
    if rms_pos == 0
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

%% features related to complex Space Vector Is :

function sv_features = extract_space_vector_features(is)


    % Extracting the metrics
    CSV_Features.Mean           = compute_mean_CSV(is);
    CSV_Features.RippleVariance = compute_envelop_ripple_variance_CSV(is);
    CSV_Features.P2PMargin      = compute_envelop_P2P_Margin_CSV(is);
    CSV_Features.Eccentricity   = compute_Geometric_Eccentricity_CSV(is);
    CSV_Features.Kurtosis       = compute_kurtosis_CSV(is);
    CSV_Features.Skewness       = compute_skewness_CSV(is);
    CSV_Features.ShapeFactor    = compute_Shape_Factor_CSV(is);
    CSV_Features.RMS            = compute_RMS_CSV(is);

    disp(struct2table(CSV_Features));
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