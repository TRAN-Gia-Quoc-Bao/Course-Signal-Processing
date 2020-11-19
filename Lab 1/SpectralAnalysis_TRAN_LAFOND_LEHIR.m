%% BE1 _ Analyse Spectrale 

% �tudiants: TRAN Gia Quoc Bao, LAFOND Arnaud, LE HIR Lena�g 
% 2e ann�e, ASI, Grenoble INP - ENSE3
% Date: 17 Avril 2020

%% Introduction: 
% L�objectif est ici de mesurer sp�cifiquement les composantes spectrales 
% des signaux afin d�avoir une id�e de la pollution harmonique produite 
% par un tel syst�me (probl�me d�analyse spectrale).

%% Commandes par d�faut
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
xlabel('�chantillons');

% On a un signal d'amplitude 15 A et une freq presque constante.

%% Calculer la valeur moyenne
moyCou = mean(crt);
% Ce signal �tant tr�s peu bruit�: moyCou = -0.078098167759585.

%% 1.1 Calculer N
% On a 65535 �chantillons et 20 cycles par 10^4 �chantillons.
% Alors 65535/(20*6.5) = 500 �chantillons/cycle
% Fe = 16384 Hz alors fondamentale 16384/500 = 32.768 Hz
% Pour choisir N en fonction de la r�solution spectrale d�sir�e,
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
ylabel('Amplitude (A�)');
xlabel('Fr�quence (Hz)');
[sp,fp] = permoy(crt, boxcar(N), 100, 20*N, fe, 'SP'); 
subplot(212);
plot(fp, sp);
grid on;
title('Permoy du signal courant avec zero-padding');
ylabel('Amplitude (A�)');
xlabel('Fr�quence (Hz)');
% Avec 0-padding de 20*N, nous avons plus de points de calcul, 
% la courbe est visuellement am�lior�e. Nous pouvons voir des 
% pics autour de la fr�quence cible qui est de 32.768 Hz.
% (il faut zoomer pour voir). Notez que l'unit� est carr�e.

%% 1.3 �chelle logarithmique
figure(3);
subplot(211);
plot(f, 10*log10(s));
grid on;
title('Permoy (log10) du signal courant sans zero-padding');
ylabel('Amplitude logarithmique (dB)');
xlabel('Fr�quence (Hz)');
subplot(212)
plot(fp, 10*log10(sp));
grid on;
title('Permoy (log10) du signal courant avec zero-padding');
ylabel('Amplitude logarithmique (dB)');
xlabel('Fr�quence (Hz)');
% Nous ne voyons que le pic � 33 Hz sur une �chelle non logarithmique.
% Le p�riodogramme en logarithme avec zero-padding est plus pr�cis, 
% notamment autour de la fr�quence fondamentale. Mais si nous ajoutons 
% trop de points, le p�riodogramme peut �tre difficile � voir. 
% Hormis le pic principal autour de 33 Hz, on voit les doubles pics autour 
% de 1200 Hz et 2400 Hz (les harmoniques) avec une modulation d'amplitude. 
% Il y a �galement 4 pics autour de 5000 Hz et une sym�trie autour de cette 
% fr�quence qui est la fr�quence de modulation.

%% 1.4 Influence des fen�tres de pond�ration
%% Fen�tre Hanning
[sh, fh] = permoy(crt, hanning(N), 100, N, fe, 'SP');
[shp, fhp] = permoy(crt, hanning(N), 100, 20*N, fe, 'SP');

figure(4);
subplot(221);
plot(fh, sh);
grid on;
title('Permoy du signal courant avec Hanning & sans zero-padding');
ylabel('Amplitude (A�)');
xlabel('Fr�quence (Hz)');
subplot(222);
plot(fhp, shp);
grid on;
title('Permoy du signal courant avec Hanning & zero-padding');
ylabel('Amplitude logarithmique (dB)');
xlabel('Fr�quence (Hz)');
subplot(223);
plot(fh, 10*log10(sh));
grid on;
title('Permoy (log10) du signal courant avec Hanning & sans zero-padding');
ylabel('Amplitude logarithmique (dB)');
xlabel('Fr�quence (Hz)');
subplot(224);
plot(fhp, 10*log10(shp));
grid on;
title('Permoy (log10) du signal courant avec Hanning & zero-padding');
ylabel('Amplitude logarithmique (dB)');
xlabel('Fr�quence (Hz)');

% Fen�tre Blackman
[sh, fh] = permoy(crt, blackman(N), 100, N, fe, 'SP');
[shp, fhp] = permoy(crt, blackman(N), 100, 20*N, fe, 'SP');

figure(6);
subplot(221);
plot(fh, sh);
grid on;
title('Permoy du signal courant avec Blackman & sans zero-padding');
ylabel('Amplitude (A�)');
xlabel('Fr�quence (Hz)');
subplot(222);
plot(fhp, shp);
grid on;
title('Permoy du signal courant avec Blackman & zero-padding');
ylabel('Amplitude logarithmique (A�)');
xlabel('Fr�quence (Hz)');
subplot(223);
plot(fh, 10*log10(sh));
grid on;
title('Permoy (log) du signal courant avec Blackman & sans zero-padding');
ylabel('Amplitude logarithmique (dB)');
xlabel('Fr�quence (Hz)');
subplot(224);
plot(fhp, 10*log10(shp));
grid on;
title('Permoy (log) du signal courant avec Blackman & zero-padding');
ylabel('Amplitude logarithmique (dB)');
xlabel('Fr�quence (Hz)');
% Il n'y a pas beaucoup de diff�rence entre les 2 fen�tres, 
% mais par rapport au cas pr�c�dent, la puissance est plus r�duite 
%(pr�s de -90 dB par rapport � -70 dB)

%% 2. Signal de vibrations statoriques
%% T�l�charger et visualiser le signal
load('signalvib.mat');

figure(7);
plot(vib);
grid on;
title('Signal vibration');
ylabel('Amplitude (m/s�)');
xlabel('�chantillons');
% On a un signal d'amplitude 0.58 m/s� et une freq presque constante

%% Calculer la valeur moyenne
moyVib = mean(vib);
% Ce signal est tr�s peu bruit�: moyVib = 1.294621486035223e-04

%% 2.1 Calculer N
% La fondamentale 12 Hz
% Pour choisir N en fonction de la r�solution spectrale d�sir�e,
% on doit avoir Fe/N << 12
% On peut prendre Fe/N = 12/10 => N = 5*Fe/6 = 10250

%% 2.2 Periodogramme

N1 = 10250;
[sv1, fv1] = permoy(vib, hanning(N1), 1, 16*N1, fe, 'SP');

figure(8);
subplot(221);
plot(fv1, sv1);
title('Permoy de vib avec Hanning & zero-padding pour N = 10250');
ylabel('Amplitude (m/s�)�');
xlabel('Fr�quence (Hz)');
subplot(222);
plot(fv1, 10*log10(sv1));
title('Permoy (log10) de vib avec Hanning & zero-padding pour N = 10250');
ylabel('Amplitude logarithmique (dB)');
xlabel('Fr�quence (Hz)');

N2 = 2500;
[sv2, fv2] = permoy(vib, hanning(N2), 1, 16*N2, fe, 'SP');
subplot(223);
plot(fv2, sv2);
title('Permoy de vib avec Hanning & zero-padding pour N = 2500');
ylabel('Amplitude (m/s�)�');
xlabel('Fr�quence (Hz)');
subplot(224);
plot(fv2, 10*log10(sv2));
title('Permoy (log10) de vib avec Hanning & zero-padding pour N = 2500');
ylabel('Amplitude logarithmique (dB)');
xlabel('Fr�quence (Hz)');

% Les variations sont trop �lev�es, alors la variance d�estimation est 
% trop grande. Il faut donc la diminuer en diminuant le d�calage et N.
% Apr�s quelques exp�riences on prend N = 2500 pour avoir un meilleur 
% compromis r�solution/variance. Il y a beaucoup moins de "variance
% d'estimation".

%% Calculer la puissance
N = 2500;
[s, f] = permoy(vib, hanning(N), 1, 16*N, fe, 'SP');

% On a essay� avec fa = 2000 Hz et fb = 6000 Hz et find mais il y a eu une erreur
% Alors on a cherch� les positions directement dans le vecteur f
% On a vu 2 kHz correspond � 6505 et 6 kHz � 19514, c'est pas mal.

puissance = sum(s(6505 : 19514)); % puissance = 0.3967 A�
pourcentage = 100*sum(s(6505 : 19514))/sum(s); % pourcentage = 47.47%
% Cela signifie que de toutes les fr�quences, celles dans la gamme 
% 2 kHz - 6 kHz repr�sentent pr�s de la moiti� de la puissance du signal.

%% Conclusion: 
% Avec les outils utilis�s tout au long des exercices, nous avons compris 
% l'utilisation de l'analyse spectrale et les liens entre la puissance 
% d'un signal et ses fr�quences.