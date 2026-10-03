clc; clear; close all;

%% 1. Charger le fichier
fprintf('Chargement de healthy_v00.mat ...\n');
data = load('faulty_stator_v00.mat');

% Vérifier la présence des variables nécessaires
if ~all(isfield(data, {'isa_f_12', 'isb_f_12', 'isc_f_12'}))
    error('Le fichier ne contient pas isa_h, isb_h, isc_h.');
end

% Extraire les courants en vecteurs lignes
isa = data.isa_f_12(:)';
isb = data.isb_f_12(:)';
isc = data.isc_f_12(:)';

fprintf('Courants chargés : isa_h, isb_h, isc_h (longueur = %d)\n', length(isa));

%% 2. Gérer le vecteur temps t
t = [];
if isfield(data, 't')
    t = data.t(:)';
    fprintf('Vecteur temps "t" trouvé (longueur %d)\n', length(t));
elseif isfield(data, 'time')
    t = data.time(:)';
    fprintf('Vecteur temps "time" trouvé (longueur %d)\n', length(t));
else
    % Générer le temps avec la fréquence d'échantillonnage
    fs = 5000;   % par défaut
    if isfield(data, 'fs') && isnumeric(data.fs) && isscalar(data.fs)
        fs = data.fs;
        fprintf('Fréquence d''échantillonnage fs = %d Hz trouvée dans le fichier\n', fs);
    else
        fprintf('fs non trouvé – utilisation de fs = %d Hz par défaut\n', fs);
    end
    N = length(isa);
    t = (0:N-1) / fs;
    fprintf('Vecteur temps t généré (longueur %d)\n', N);
end

%% 3. Supprimer les 2 premières secondes (garder les 8 dernières)
% Trouver l'index correspondant à t >= 2 secondes
idx = (t >= 2);   % indices logiques à partir de 2s

% Appliquer la coupure à tous les signaux et au temps
t = t(idx);
isa_crop = isa(idx);
isb_crop = isb(idx);
isc_crop = isc(idx);

fprintf('Suppression des 2 premières secondes : %d échantillons conservés (sur %d)\n', ...
        length(t), length(t));

%% 4. Calculer le vecteur espace i_s à partir des signaux coupés
K = 2/3;   % gain puissance-invariante (mettre 1 pour amplitude-invariante)
alpha = exp(1j * 2*pi/3);
beta  = exp(1j * 4*pi/3);    % identique à alpha^2

i_s = K * (isa_crop + alpha * isb_crop + beta * isc_crop);

fprintf('i_s calculé (vecteur complexe, longueur %d)\n', length(i_s));

%% 5. Sauvegarder dans signals_data.mat (uniquement i_s et t_crop)
signals = struct();
signals.i_s = i_s;

output_filename = 'signals_data.mat';
save(output_filename, 'signals', 't');
fprintf('\nSauvegardé dans "%s" :\n', output_filename);
fprintf('  - signals.i_s : vecteur complexe [1 x %d]\n', length(i_s));
fprintf('  - t_crop      : vecteur temps [1 x %d] (8 secondes)\n', length(t));

%% 6. Tracer le signal complexe i_s (parties réelle et imaginaire)
figure('Name', 'Vecteur espace courant i_s (8 dernières secondes)', ...
       'Position', [100 100 800 600]);

subplot(2,1,1);
plot(t, real(i_s), 'b', 'LineWidth', 1.2);
xlabel('Temps (s)');
ylabel('Re(i_s)');
title('Partie réelle de \vec{i}_s');
grid on;
xlim([t(1) t(end)]);

subplot(2,1,2);
plot(t, imag(i_s), 'r', 'LineWidth', 1.2);
xlabel('Temps (s)');
ylabel('Im(i_s)');
title('Partie imaginaire de \vec{i}_s');
grid on;
xlim([t(1) t(end)]);


fprintf('\nTracé terminé.\n');