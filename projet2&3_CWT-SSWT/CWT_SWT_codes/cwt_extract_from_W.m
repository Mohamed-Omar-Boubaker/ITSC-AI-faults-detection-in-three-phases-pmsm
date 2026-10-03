function [x_rec, band_info] = cwt_extract_from_W(fe_k, p, W, f)
    page = 1 + (fe_k <= 0);
    inverse_page = 3 - page;  % 1->2, 2->1
    [~, idx_start] = min(abs(f - p.userBand(1)));
    [~, idx_end]   = min(abs(f - p.userBand(2)));
    idx = min(idx_start, idx_end) : max(idx_start, idx_end);

    W_filt = zeros(size(W));
    W_filt(:, :, inverse_page) = 0;
    mask = zeros(size(W,1), size(W,2));
    mask(idx, :) = 1;
    W_filt(:, :, page) = W(:, :, page) .* mask;

    f_min = min(f(idx));
    f_max = max(f(idx));

    x_rec = icwt(W_filt, p.wavelet, ...
                 'WaveletParameters', [p.gamma, p.gamma*p.beta_opt], ...
                 'VoicesPerOctave', p.voices_per_octave);

    band_info.idx    = idx;
    band_info.f_low  = f_min;
    band_info.f_high = f_max;
    band_info.page   = page;
end
