%% generate_itsc_signal.m
% Generates synthetic three-phase stator currents of an induction motor with an
% inter-turn short circuit (ITSC) in phase A and saves them to a .mat file.
% Saved variables: ia, ib, ic, fs  (+ severity_percent)  -> matches the default
% cfg.VarNames = {'ia','ib','ic'} of MatSignalLoader, so no config change needed.

clear; clc;

%% ===== USER SETTINGS =====
severity_percent = 3;      % ITSC severity: 0 .. 30 (% of shorted turns in phase A)
fs        = 10000;          % sampling rate [Hz]
duration  = 10;              % signal length [s]
f_supply  = 200;             % supply frequency [Hz]
I_rms     = 0.85;             % healthy rms phase current [A]
slip      = 0.03;           % rotor slip
pole_pairs = 2;             % pole pairs
noise_snr_dB = 40;          % measurement noise (set Inf for none)
out_dir  = 'C:\Users\moham\OneDrive\Desktop\SeverityPipeline';
out_file = fullfile(out_dir, sprintf('itsc_%02dpct.mat', round(severity_percent)));%% =========================

assert(severity_percent >= 0 && severity_percent <= 30, 'Severity must be within 0..30 %%.');

sev = severity_percent / 100;           % 0 .. 0.30
t   = (0:1/fs:duration - 1/fs).';
w   = 2*pi*f_supply;
Ipk = sqrt(2) * I_rms;
shift = [0, -2*pi/3, 2*pi/3];           % a, b, c

%% Healthy positive-sequence currents
I = zeros(numel(t), 3);
for k = 1:3
    I(:, k) = Ipk * cos(w*t + shift(k));
end

%% ITSC effects
% 1) Negative-sequence component (main unbalance signature), aligned so phase A grows
k_neg = 1.2;                            % gain: negative-seq amplitude / Ipk per unit severity
In    = k_neg * sev * Ipk;
for k = 1:3
    I(:, k) = I(:, k) + In * cos(w*t - shift(k));   % reversed phase sequence
end

% 2) Extra current in the faulted phase (shorted-turn circulating effect)
I(:, 1) = I(:, 1) + 0.5 * sev * Ipk * cos(w*t + deg2rad(20));

% 3) Fault harmonics: f_sc = f_s * ( n*(1-s)/p  +/-  k ),  n = 1,3 ; k = 1,3
amp0 = 0.25 * sev * Ipk;                % harmonic amplitude scales with severity
n_list = [1 3];  k_list = [1 3];
for n = n_list
    for kk = k_list
        for sgn = [-1 1]
            f_sc = f_supply * abs(n*(1 - slip)/pole_pairs + sgn*kk);
            a    = amp0 / (n*kk);
            ph   = 2*pi*rand;
            for ph_idx = 1:3
                % faulted phase gets the strongest component
                g = 1; if ph_idx > 1, g = 0.6; end
                I(:, ph_idx) = I(:, ph_idx) + g * a * cos(2*pi*f_sc*t + ph + shift(ph_idx));
            end
        end
    end
end

%% Measurement noise
if isfinite(noise_snr_dB)
    for k = 1:3
        p_sig = mean(I(:, k).^2);
        I(:, k) = I(:, k) + sqrt(p_sig / 10^(noise_snr_dB/10)) * randn(size(t));
    end
end

ia = I(:, 1);  ib = I(:, 2);  ic = I(:, 3);

%% Save
save(out_file, 'ia', 'ib', 'ic', 'fs', 'severity_percent');
fprintf('Saved %s  (severity = %.1f %%, fs = %d Hz, %d samples)\n', ...
        out_file, severity_percent, fs, numel(t));

%% Quick look
figure('Name', sprintf('ITSC %.0f%%', severity_percent));
nShow = round(0.06 * fs);
subplot(2,1,1); plot(t(1:nShow), I(1:nShow, :)); grid on;
legend('ia','ib','ic'); xlabel('t [s]'); ylabel('A'); title('Phase currents (first 60 ms)');
subplot(2,1,2);
N = numel(t); F = abs(fft(ia .* hann(N))) / N * 4; f = (0:N-1).' * fs / N;
semilogy(f(1:N/2), F(1:N/2) + eps); xlim([0 500]); grid on;
xlabel('Hz'); ylabel('|Ia|'); title('Spectrum of ia');
