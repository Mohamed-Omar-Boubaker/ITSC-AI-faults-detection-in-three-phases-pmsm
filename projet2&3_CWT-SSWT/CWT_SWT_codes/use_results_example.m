% use_results_example.m
% Example of loading and using pipeline_results.mat in a separate script.
% Run Wavelet_Spectrum_Shifter() first to generate the file.

clear; clc;

% ---- 1. Load ----
loaded = load('pipeline_results.mat');   % returns a struct with field 'results'
results = loaded.results;

% ---- 2. List datasets available ----
sig_names = fieldnames(results);
disp('Datasets found:');
disp(sig_names);

% ---- 3. Pick one dataset ----
sig_name = sig_names{1};
R = results.(sig_name);

% ---- 4. Access top-level fields ----
x_original = R.original_signal;
x_rec_total = R.x_rec_total;
Fs = R.Fs;

% ---- 5. Loop over each detected-frequency entry ----
for k = 1:numel(R.cwt_info)
    entry = R.cwt_info{k};

    fprintf('Freq = %+.2f Hz | amp = %.4f | beta = %.4f | vpo = %d\n', ...
            entry.freq, entry.amp, entry.beta, entry.vpo);

    % Scalogram and frequency axis
    W = entry.W;
    f = entry.f;

    % Bands actually used for reconstruction
    band_pos = entry.band_pos;
    band_neg = entry.band_neg;

    % Reconstructed components
    x_rec_pos = entry.x_rec_pos;
    x_rec_neg = entry.x_rec_neg;

    % ... do whatever downstream processing you need here ...
end
