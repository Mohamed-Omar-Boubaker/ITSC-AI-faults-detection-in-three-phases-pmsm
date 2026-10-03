function feats = extract_cwt_features(results)

    W = results.cwt_info.W;
    f = results.cwt_info.f;
    
    feats.freq              = entry.freq;
    feats.total_energy      = compute_total_Scalogram_Energy(W);
    feats.positive_energy   = compute_total_Positive_Scalogram_Energy(W);
    feats.negative_energy   = compute_total_Negative_Scalogram_Energy(W);
    feats.energy_ratio      = compute_Scalogram_Energy_Ratio(W);
    feats.positive_entropy  = compute_Positive_Shannon_Entropy(W);
    feats.negative_entropy  = compute_Negative_Shannon_Entropy(W);
    feats.spectral_centroid = compute_Spectral_Center_of_Gravity(W, f);
   
end

%% ── LOCAL METRIC FUNCTIONS ─────────────────────────────────────────
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