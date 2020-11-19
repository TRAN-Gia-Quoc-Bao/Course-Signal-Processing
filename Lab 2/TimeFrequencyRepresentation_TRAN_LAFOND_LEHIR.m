%% BE2 _ Représentations Temps-Fréquence

% Étudiants: TRAN Gia Quoc Bao, LAFOND Arnaud, LE HIR Lenaïg 
% 2e année, ASI, Grenoble INP - ENSE3
% Date: 17 Avril 2020

%% Introduction: 
% Dans ce BE, le spectrogramme sera utilisé pour les représentations 
% temps-fréquence pour l'analyse spectrale des signaux non stationnaires.

%% Commandes par défault
close all;
clear all;
clc;

%% Telecharger et visualiser le signal
load('signalcrtns.mat');

figure(1);
plot(crtns);
grid on;
title('Signal de courant de phase de la machine asynchrone');
xlabel('Échantillons');
ylabel('Amplitude (A)');

% On voit 4 phases (en échantillons):
% 1 - environ 23000: freq = 2*475/60 = 15.8 Hz, amplitude const = 2 A
% 23000 - 65000: freq augmente, amplitude const = 2.4 A
% 65000 - 107000: freq diminue, amplitude diminue de 2 à 1.8 A
% reste: comme phase 1

%% 1. Analyse autour du fondamental 

% On determine N :
% 1er phase : 475 echantillon dans une periode
% donc fréquence du courant electrique dans la premiere phase est de
% 475*2/60 = 15.8 Hz
% d'apres la courbe de consigne : delta_consigne/delta_t = 1000-475 = 525tr/min
% or delta_f = delta_consigne*2/60
% donc delta_f/delta_t = 525*2/60 = 17.5 Hz

Nopt = fe*sqrt(1/17.5); % Nopt = 3.9165e+03
N = 4000; % on prend cette valeur

[stft, f, t] = spectrogram(crtns, hanning(N), ceil(0.5*N), 8*N, fe);
rtf = abs(stft).^2;

figure(2);
imagesc(t, f, rtf); 
grid on;
axis xy; 
colorbar;
xlim([0 7.3]);
ylim([0 100]);
title('Spectrogramme');

figure(3);
imagesc(t, f, 10*log10(rtf)); 
grid on; 
axis xy; 
colorbar;
%zoom
%xlim([0 7.3]);
%ylim([0 100]);
title('Spectrogramme échelle Log');

% Nous avons essayé avec différentes fenêtres et nous pensons que 
% Hanning donne les résultats les plus clairs.

% La reponse est un peu differente de la commande car le systeme répond
% lentement. La vitesse de rotation suit la consigne impose mais avec un
% retard en fait il repond comme un systeme du premiere ordre.
% Pour ameliorer, il faudrait modofier en ordre superieur.

% Le zero padding ne changepas grand chose au resultat.

%% 2. Analyse sur toute la bande fréquentielle disponible
% Autour de 0 Hz, on retrouve des raies que nous supposons correspondre à  
% la fréquence de la phase 1 (16 Hz).

% Autour de 500-600Hz, on peut penser au phénomène d'encoche avec 
% les harmoniques situé vers 1200 puis 1800 Hz.

% Vers 5000Hz, on retrouve les raies correspondante à la commande MLI avec
% la aussi le phénomène d'encoche.

%% 3. Analyse autour du phénomène d'encoches 

% On se situe au tour du 72*vr, alors on doit recalculer N
N = floor(Nopt/sqrt(72));

% Nous avons essayé avec différentes fenêtres et nous pensons que 
% Hanning donne les résultats les plus clairs.

[stft, f, t] = spectrogram(crtns, hanning(N), floor(0.7*N), N, fe);
rtf = abs(stft).^2;

figure(4);
imagesc(t, f, rtf);
grid on;
axis('xy');
colorbar;
xlim([0 7.3]);
ylim([0 100]);
title('Spectrogramme');

figure(5);
imagesc(t, f, 10*log10(rtf));
axis('xy');
title('Analyse spectrale temps-fréquence autour du phénomène d encoches du signal');

% Seule la clarté est meilleure car nous avons optimisé N pour une 
% détection du phénomène des crans à 550-600 Hz avec ses harmoniques et 
% autour de la fréquence du MLI.

%% Conclusion: 
% Le spectrogramme, avec tous les outils vus dans le BE, est utile pour 
% l'analyse spectrale de signaux non stationnaires. Pour l'utiliser, 
% d'abord une analyse des différentes phases du signal, puis de toute 
% la bande fréquentielle disponible, et on peut avoir la bonne zone de fréquences.