function print_summary(results)
    fprintf('\n\n==================== FINAL SUMMARY ====================\n');

    fprintf('Original signal: length = %d, Fs = %.1f Hz\n', ...
            length(results.original_signal), results.Fs);

    det = results.detected_freqs;
    fprintf('Detected FFT peaks (freq, amp):\n');
    for i = 1:size(det,1)
        fprintf('  %+.2f Hz  ->  %.4f\n', det(i,1), det(i,2));
    end

    fprintf('\nCWT entries (Targeted Peak Components):\n');
    
    % Accessing directly as a flat struct field (no loop, no cell array brackets)
    entry = results.cwt_info;
    
    % FIXED: Updated field names to match your new freq_pos, freq_neg, amp_pos, amp_neg fields
    fprintf('  Targeted Pos Freq (fe_pos) = %+.2f Hz, amp = %.4f\n', entry.freq_pos, entry.amp_pos);
    fprintf('  Targeted Neg Freq (fe_neg) = %+.2f Hz, amp = %.4f\n', entry.freq_neg, entry.amp_neg);
    fprintf('  Wavelet Parameters: beta = %.4f, VPO = %d\n', entry.beta, entry.vpo);
    
    fprintf('    Band positive: [%.2f, %.2f] Hz (page %d)\n', ...
            entry.band_pos.f_low, entry.band_pos.f_high, entry.band_pos.page);
    fprintf('    Band negative: [%.2f, %.2f] Hz (page %d)\n', ...
            entry.band_neg.f_low, entry.band_neg.f_high, entry.band_neg.page);
            
    fprintf('    Scalogram size: %d x %d x %d\n', size(entry.W,1), size(entry.W,2), size(entry.W,3));
    fprintf('    Reconstructed signals: pos length %d, neg length %d\n', ...
            length(entry.x_rec_pos), length(entry.x_rec_neg));
            
    fprintf('    Total Reconstructed Signal Length: %d\n', length(results.x_rec_total));

    fprintf('========================================================\n');
end