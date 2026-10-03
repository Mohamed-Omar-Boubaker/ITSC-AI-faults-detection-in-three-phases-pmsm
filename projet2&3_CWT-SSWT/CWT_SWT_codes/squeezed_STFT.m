%% Synchrosqueezing Transform (SST) fait maison — sans Wavelet Toolbox
% Basee sur la STFT. Fonctionne pour signaux reels ET complexes.
% Reference algorithme : Daubechies, Lu, Wu (2011) ; Thakur & Wu (2013)
data = load('faulty_stator_v00.mat'); 

N=length(data.Te_f_12);
% keep only the last 8s 
fs=5000;
t = (0:N-1) / fs;
idx = (t >= 2);   
t = t(idx);
    

    
% Extraction des champs requis
Te = data.Te_f_12(:, 1);
ia = data.isa_f_12(:, 1);
ib = data.isb_f_12(:, 1);
ic = data.isc_f_12(:, 1);
    
% keep only useful 8s
Te = Te(idx);
ia = ia(idx);
ib = ib(idx);
ic = ic(idx);
     
% Compute the complex space vector 
x = sqrt(2/3) * (ia + ib * exp(1i*2*pi/3) + ic * exp(1i*4*pi/3));

N=length(x);


% %% 1. Construction d'un signal complexe (chirp I/Q)
% fs = 1000;              % frequence d'echantillonnage (Hz)
% t  = (0:1/fs:2-1/fs)';  % duree 2 s
% 
% % Chirp complexe : frequence instantanee variant de 50 Hz a 200 Hz
% f0 = 50; f1 = 200; T = t(end);
% k  = (f1 - f0) / T;                     % taux de balayage (Hz/s)
% phase = pi*(f0*t + k/2*t.^2);
% x1 = exp(1i*phase);                     % composante 1 : chirp montant
% % x1 = 5*exp(1i*2*pi*50*t);
% % Deuxieme composante : frequence NEGATIVE constante (impossible a
% % distinguer avec un signal reel, mais naturel avec un signal complexe)
% x2 = 0.6*exp(-1i*2*pi*80*t);
% 
% x = x1;                            % signal complexe final (bruit ajoute)
% %x = x + 0.02*(randn(size(x)) + 1i*randn(size(x)));

%% 2. STFT complexe (fenetre glissante + FFT complete, freq. + et -)
win_len = 256;
hop     = 1;
nfft    = 512;    
win     = hamming(win_len);

n_frames = floor((length(x) - win_len)/hop) + 1;
STFT     = zeros(nfft, n_frames);
dSTFT    = zeros(nfft, n_frames);   % derivee temporelle de la STFT

% frequences FFT centrees (positives ET negatives -> essentiel en complexe)
freqs = (-nfft/2:nfft/2-1)' * (fs/nfft);

for m = 1:n_frames
    idx = (m-1)*hop + (1:win_len);
    seg = x(idx) .* win;

    X  = fftshift(fft(seg, nfft));
    % derivee de la fenetre pour estimer d/dt du signal fenetre
    dwin = gradient(win) * fs;
    dseg = x(idx) .* dwin;
    dX = fftshift(fft(dseg, nfft));
    
    
    
    STFT(:,m)  = X;
    dSTFT(:,m) = dX;
end

time_frames = ((0:n_frames-1)*hop + win_len/2) / fs;

%% 3. Estimation de la frequence instantanee et reassignation (squeezing)
threshold = 1e-8 * max(abs(STFT(:)));   % seuil pour eviter division par ~0
omega = nan(size(STFT));
mask = abs(STFT) > threshold;

% IMPORTANT : il faut ajouter la frequence du bin d'analyse f !
% Formule correcte (Auger-Flandrin / Daubechies-Maes) :
%   omega_hat(t,f) = f - (1/2pi) * Im[ S_h'(t,f) / S_h(t,f) ]
freq_grid = repmat(freqs, 1, n_frames);   % meme taille que STFT
omega(mask) = freq_grid(mask) - imag(dSTFT(mask) ./ STFT(mask)) / (2*pi);  % Hz

% Grille de sortie synchrosqueezee (memes bornes de frequence)
SST = zeros(size(STFT));
df  = freqs(2) - freqs(1);
a=0;
for m = 1:n_frames
    for k_bin = 1:nfft
        if mask(k_bin, m)
            f_est = omega(k_bin, m);
            % on ignore les frequences estimees hors de la bande utile
            if f_est >= freqs(1) && f_est <= freqs(end)
                [~, k_target] = min(abs(freqs - f_est));
                SST(k_target, m) = SST(k_target, m) + STFT(k_bin, m) * df;
            end
        end
        a=a+1;
    end
end

%% 4. Affichage : STFT classique vs SST — échelle de couleur commune
% Calcul des amplitudes
amp_STFT = abs(STFT)/sum(win);
amp_SST   = abs(SST)/sum(win);

% Déterminer la valeur maximale commune (min = 0)
max_common = max([amp_STFT(:); amp_SST(:)]);

figure();

subplot(2,1,1);
imagesc(time_frames, freqs, amp_STFT);
axis xy;
colormap(jet);
caxis([0 max_common]);   % ← Échelle unifiée
colorbar;
xlabel('Temps (s)'); ylabel('Fréquence (Hz)');
title('STFT classique (signal complexe : freq. + et - visibles)');

subplot(2,1,2);
imagesc(time_frames, freqs, amp_SST);
axis xy;
colormap(jet);
caxis([0 max_common]);   % ← Même échelle
colorbar;
xlabel('Temps (s)'); ylabel('Fréquence (Hz)');
title('Synchrosqueezing Transform (SST) — résolution améliorée');

%% 5. Verification : extraction d'une ridge (optionnel)
% Exemple simple : recherche du maximum d'energie a chaque instant
[~, ridge_idx] = max(abs(SST), [], 1);
ridge_freq = freqs(ridge_idx);

figure;
plot(time_frames, ridge_freq, 'LineWidth', 1.5);
xlabel('Temps (s)'); ylabel('Frequence dominante estimee (Hz)');
title('Ridge extraite de la SST (composante d''energie max)');
grid on;


%% 6. Reconstruction du signal à partir de la STFT classique
x_rec_STFT = zeros(size(x));
win_sum_sq = zeros(size(x));   % accumulateur pour la normalisation

for m = 1:n_frames
    % Segment temporel correspondant
    idx = (m-1)*hop + (1:win_len);
    
    % Colonne STFT -> domaine temporel
    seg = ifft(ifftshift(STFT(:, m)));   % longueur nfft, complexe
    seg = seg(1:win_len);                % on ne garde que la partie non nulle
    
    % Synthèse par fenêtrage et recouvrement (overlap-add)
    x_rec_STFT(idx) = x_rec_STFT(idx) + seg .* win;
    win_sum_sq(idx)  = win_sum_sq(idx)  + win.^2;
end

% Normalisation (compense le chevauchement)
x_rec_STFT = x_rec_STFT ./ max(win_sum_sq, eps);

err_STFT = norm(x - x_rec_STFT)^2 / norm(x)^2;
fprintf('Erreur relative de reconstruction STFT : %.2e\n', err_STFT);

%% 7. Reconstruction du signal à partir de la SST (corrigé)
x_rec_SST = zeros(size(x));
win_sum_sq = zeros(size(x));   % normalisation par fenêtre²

for m = 1:n_frames
    idx = (m-1)*hop + (1:win_len);
    
    % Retirer le facteur df pour revenir aux coefficients de Fourier
    seg = ifft(ifftshift(SST(:, m)/df ));
    seg = seg(1:win_len);   % on jette le zéro‑padding
    
    % Overlap‑add identique à la STFT
    x_rec_SST(idx) = x_rec_SST(idx) + seg .* win;
    win_sum_sq(idx) = win_sum_sq(idx) + win.^2;
end

x_rec_SST = x_rec_SST ./ max(win_sum_sq, eps);

err_SST = norm(x - x_rec_SST)^2 / norm(x)^2;
fprintf('Erreur relative de reconstruction SST  : %.2e\n', err_SST);

%% 8. Visualisation des reconstructions (optionnel)
figure;
subplot(3,1,1);
plot(t, real(x), 'b', t, imag(x), 'r--');
title('Signal original (réel et imaginaire)'); xlabel('Temps (s)'); legend('Réel','Imag');
xlim([2.4 2.5]);   % zoom sur le début pour voir les détails

subplot(3,1,2);
plot(t, real(x_rec_STFT), 'b', t, imag(x_rec_STFT), 'r--');
title(sprintf('Reconstruction depuis STFT (erreur rel = %.2e)', err_STFT));
xlabel('Temps (s)'); xlim([2.4 2.5]);

subplot(3,1,3);
plot(t, real(x_rec_SST), 'b', t, imag(x_rec_SST), 'r--');
title(sprintf('Reconstruction depuis SST (erreur rel = %.2e)', err_SST));
xlabel('Temps (s)'); xlim([2.4 2.5]);


%% 
% 1. Conservation de l'énergie réelle / Module (Norme L1)
E_STFT_module = sum(abs(STFT(:))) * df;
E_SST_module  = sum(abs(SST(:))); 

% 2. Conservation de l'énergie quadratique / Puissance (Norme L2)
E_STFT_energie = sum(abs(STFT(:)).^2) * df;
E_SST_energie  = sum(abs(SST(:)).^2);

disp(['Énergie STFT : ', num2str(E_STFT_energie)]);
disp(['Énergie SST  : ', num2str(E_SST_energie)]);

disp(['module STFT : ', num2str(E_STFT_module)]);
disp(['module SST  : ', num2str(E_SST_module)]);


