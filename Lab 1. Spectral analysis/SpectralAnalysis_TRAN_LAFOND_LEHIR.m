%% BE1 _ Analyse Spectrale 

% Étudiants: TRAN Gia Quoc Bao, LAFOND Arnaud, LE HIR Lenaïg 
% 2e année, ASI, Grenoble INP - ENSE3
% Date: 17 Avril 2020

%% Introduction: 
% L’objectif est ici de mesurer spécifiquement les composantes spectrales 
% des signaux afin d’avoir une idée de la pollution harmonique produite 
% par un tel système (problème d’analyse spectrale).

%% Commandes par défaut
close all;
clear all;
clc;

%% 1. Signal the courant statorique
%% Telecharger et visualiser le signal
load('signalcrt.mat');

figure(1);
plot(crt);
grid on;
title('Signal courant');
ylabel('Amplitude (A)');
xlabel('Échantillons');

% On a un signal d'amplitude 15 A et une freq presque constante.

%% Calculer la valeur moyenne
moyCou = mean(crt);
% Ce signal étant très peu bruité: moyCou = -0.078098167759585.

%% 1.1 Calculer N
% On a 65535 échantillons et 20 cycles par 10^4 échantillons.
% Alors 65535/(20*6.5) = 500 échantillons/cycle
% Fe = 16384 Hz alors fondamentale 16384/500 = 32.768 Hz
% Pour choisir N en fonction de la résolution spectrale désirée,
% on doit avoir Fe/N << 32.768
% On peut prendre Fe/N = 32.768/10 => N = Fe/3 environ 5500

%% 1.2 Periodogramme
N = 5500;

figure(2);
[s, f] = permoy(crt, boxcar(N), 100, N, fe, 'SP'); 
subplot(211);
plot(f, s);
grid on;
title('Permoy du signal courant sans zero-padding');
ylabel('Amplitude (A²)');
xlabel('Fréquence (Hz)');
[sp,fp] = permoy(crt, boxcar(N), 100, 20*N, fe, 'SP'); 
subplot(212);
plot(fp, sp);
grid on;
title('Permoy du signal courant avec zero-padding');
ylabel('Amplitude (A²)');
xlabel('Fréquence (Hz)');
% Avec 0-padding de 20*N, nous avons plus de points de calcul, 
% la courbe est visuellement améliorée. Nous pouvons voir des 
% pics autour de la fréquence cible qui est de 32.768 Hz.
% (il faut zoomer pour voir). Notez que l'unité est carrée.

%% 1.3 Échelle logarithmique
figure(3);
subplot(211);
plot(f, 10*log10(s));
grid on;
title('Permoy (log10) du signal courant sans zero-padding');
ylabel('Amplitude logarithmique (dB)');
xlabel('Fréquence (Hz)');
subplot(212)
plot(fp, 10*log10(sp));
grid on;
title('Permoy (log10) du signal courant avec zero-padding');
ylabel('Amplitude logarithmique (dB)');
xlabel('Fréquence (Hz)');
% Nous ne voyons que le pic à 33 Hz sur une échelle non logarithmique.
% Le périodogramme en logarithme avec zero-padding est plus précis, 
% notamment autour de la fréquence fondamentale. Mais si nous ajoutons 
% trop de points, le périodogramme peut être difficile à voir. 
% Hormis le pic principal autour de 33 Hz, on voit les doubles pics autour 
% de 1200 Hz et 2400 Hz (les harmoniques) avec une modulation d'amplitude. 
% Il y a également 4 pics autour de 5000 Hz et une symétrie autour de cette 
% fréquence qui est la fréquence de modulation.

%% 1.4 Influence des fenêtres de pondération
%% Fenêtre Hanning
[sh, fh] = permoy(crt, hanning(N), 100, N, fe, 'SP');
[shp, fhp] = permoy(crt, hanning(N), 100, 20*N, fe, 'SP');

figure(4);
subplot(221);
plot(fh, sh);
grid on;
title('Permoy du signal courant avec Hanning & sans zero-padding');
ylabel('Amplitude (A²)');
xlabel('Fréquence (Hz)');
subplot(222);
plot(fhp, shp);
grid on;
title('Permoy du signal courant avec Hanning & zero-padding');
ylabel('Amplitude logarithmique (dB)');
xlabel('Fréquence (Hz)');
subplot(223);
plot(fh, 10*log10(sh));
grid on;
title('Permoy (log10) du signal courant avec Hanning & sans zero-padding');
ylabel('Amplitude logarithmique (dB)');
xlabel('Fréquence (Hz)');
subplot(224);
plot(fhp, 10*log10(shp));
grid on;
title('Permoy (log10) du signal courant avec Hanning & zero-padding');
ylabel('Amplitude logarithmique (dB)');
xlabel('Fréquence (Hz)');

% Fenêtre Blackman
[sh, fh] = permoy(crt, blackman(N), 100, N, fe, 'SP');
[shp, fhp] = permoy(crt, blackman(N), 100, 20*N, fe, 'SP');

figure(6);
subplot(221);
plot(fh, sh);
grid on;
title('Permoy du signal courant avec Blackman & sans zero-padding');
ylabel('Amplitude (A²)');
xlabel('Fréquence (Hz)');
subplot(222);
plot(fhp, shp);
grid on;
title('Permoy du signal courant avec Blackman & zero-padding');
ylabel('Amplitude logarithmique (A²)');
xlabel('Fréquence (Hz)');
subplot(223);
plot(fh, 10*log10(sh));
grid on;
title('Permoy (log) du signal courant avec Blackman & sans zero-padding');
ylabel('Amplitude logarithmique (dB)');
xlabel('Fréquence (Hz)');
subplot(224);
plot(fhp, 10*log10(shp));
grid on;
title('Permoy (log) du signal courant avec Blackman & zero-padding');
ylabel('Amplitude logarithmique (dB)');
xlabel('Fréquence (Hz)');
% Il n'y a pas beaucoup de différence entre les 2 fenêtres, 
% mais par rapport au cas précédent, la puissance est plus réduite 
%(près de -90 dB par rapport à -70 dB)

%% 2. Signal de vibrations statoriques
%% Télécharger et visualiser le signal
load('signalvib.mat');

figure(7);
plot(vib);
grid on;
title('Signal vibration');
ylabel('Amplitude (m/s²)');
xlabel('Échantillons');
% On a un signal d'amplitude 0.58 m/s² et une freq presque constante

%% Calculer la valeur moyenne
moyVib = mean(vib);
% Ce signal est très peu bruité: moyVib = 1.294621486035223e-04

%% 2.1 Calculer N
% La fondamentale 12 Hz
% Pour choisir N en fonction de la résolution spectrale désirée,
% on doit avoir Fe/N << 12
% On peut prendre Fe/N = 12/10 => N = 5*Fe/6 = 10250

%% 2.2 Periodogramme

N1 = 10250;
[sv1, fv1] = permoy(vib, hanning(N1), 1, 16*N1, fe, 'SP');

figure(8);
subplot(221);
plot(fv1, sv1);
title('Permoy de vib avec Hanning & zero-padding pour N = 10250');
ylabel('Amplitude (m/s²)²');
xlabel('Fréquence (Hz)');
subplot(222);
plot(fv1, 10*log10(sv1));
title('Permoy (log10) de vib avec Hanning & zero-padding pour N = 10250');
ylabel('Amplitude logarithmique (dB)');
xlabel('Fréquence (Hz)');

N2 = 2500;
[sv2, fv2] = permoy(vib, hanning(N2), 1, 16*N2, fe, 'SP');
subplot(223);
plot(fv2, sv2);
title('Permoy de vib avec Hanning & zero-padding pour N = 2500');
ylabel('Amplitude (m/s²)²');
xlabel('Fréquence (Hz)');
subplot(224);
plot(fv2, 10*log10(sv2));
title('Permoy (log10) de vib avec Hanning & zero-padding pour N = 2500');
ylabel('Amplitude logarithmique (dB)');
xlabel('Fréquence (Hz)');

% Les variations sont trop élevées, alors la variance d’estimation est 
% trop grande. Il faut donc la diminuer en diminuant le décalage et N.
% Après quelques expériences on prend N = 2500 pour avoir un meilleur 
% compromis résolution/variance. Il y a beaucoup moins de "variance
% d'estimation".

%% Calculer la puissance
N = 2500;
[s, f] = permoy(vib, hanning(N), 1, 16*N, fe, 'SP');

% On a essayé avec fa = 2000 Hz et fb = 6000 Hz et find mais il y a eu une erreur
% Alors on a cherché les positions directement dans le vecteur f
% On a vu 2 kHz correspond à 6505 et 6 kHz à 19514, c'est pas mal.

puissance = sum(s(6505 : 19514)); % puissance = 0.3967 A²
pourcentage = 100*sum(s(6505 : 19514))/sum(s); % pourcentage = 47.47%
% Cela signifie que de toutes les fréquences, celles dans la gamme 
% 2 kHz - 6 kHz représentent près de la moitié de la puissance du signal.

%% Conclusion: 
% Avec les outils utilisés tout au long des exercices, nous avons compris 
% l'utilisation de l'analyse spectrale et les liens entre la puissance 
% d'un signal et ses fréquences.