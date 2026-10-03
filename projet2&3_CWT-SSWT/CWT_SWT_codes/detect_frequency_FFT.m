function [freqs, amps, fax, fft_orig] = detect_frequency_FFT(sig_x, p)
    N = length(sig_x);
    fft_orig = abs(fftshift(fft(sig_x))/N);
    fax = (-N/2 : N/2-1) * (p.fs/N);
    win = abs(fax) <= p.f_disp;
    f_search = fax(win);
    mag_search = fft_orig(win);
    [peaks, locs] = findpeaks(mag_search, f_search, ...
        'MinPeakHeight', p.peak_thr * max(fft_orig), ...
        'MinPeakDistance', p.peak_min_dist);
    freqs = locs(:)';
    amps = peaks(:)';
end
