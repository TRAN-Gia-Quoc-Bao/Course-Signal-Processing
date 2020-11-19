%% BE2 _ Optimal Linear Filtering _ TRAN Gia Quoc Bao, ASI 2nd year

%% I. Introduction
% In this BE we will learn how to get the fetus' cardio signal out of the 
% combination fetus + mother cardio signal by applying the Widrow filter
% and the "noise only" signal (the mother's). 

%% Default commands
close all;
clear all;
clc;

%% II. The Widrow experiment:
% The major difficulty when estimating the fetal ECG lies in the fact that 
% its amplitude is 2-1000 times lower than maternal ECG (smaller heart / 
% environments to cross). Direct processing of records is therefore difficult; 
% it is necessary to eliminate interference from the maternal ECG. We'll do that
% using the Widrow algorithm (LMS filter).

%% Data loading

% The description :
% The data include the ECG signals recorded on a pregnant woman (8 electrodes)
% Sampling: Fe = 250 Hz
% Nb of samples: 2500, i.e. an observation time of 10s
% Formatting:
% Each line of the 8 x 2500 fetal_ecg matrix corresponds to
% one of the 8 electrodes:
% - electrodes 1-5: abdomen, ECG of the mother + fetus -> channel signal
% - electrodes 6-8: chest, only the mother's heart rate is measurement -> "noise only" references

load foetal_ecg.mat;   % load data
Fs = 250; % sampling frequency in Hz
N = size(foetal_ecg, 2); % number of samples
time = (0 : N - 1)/Fs; % discrete time

%% Signal shaping

% Signal to denoise: abdomen electrode number 'elec_abdo'
% Rq: the signals are already centered

elec_abdo = 1;
abdomen_ecg = foetal_ecg(elec_abdo, :); % line vector 1 x N

% "Noise only" references (secondary channel)
% Remark: the signals are already centered

ireferences = [6]; % the reference electrode (s) / signals
nRef = length(ireferences) ; % number of reference signals
thorax_ecg = foetal_ecg(ireferences, 1 : N) ; % nRef x N <- multi-ref

clear foetal_ecg ;

%% Visualization of signals

figure('Name', 'Visualization of waveforms', 'NumberTitle', 'off');
subplot(211);
plot(time, thorax_ecg'); 
grid on ; 
xlim([0 10]);
xlabel('Time (seconds)');
ylabel('Signal (voltage)');
title(['Reference(s) : Thoracic electrode number ', int2str(ireferences)]);
subplot(212);
plot(time, abdomen_ecg, '-b') ; 
grid on ; 
xlim([0 10]);
xlabel('Time (seconds)') ;
ylabel('Signal (voltage)');
title(['Signal : Abdomen electrode number ' int2str(elec_abdo)]);
% From here we see that the abdomen contains a lot of background noises
% while for the thoracic it's almost purely the mother's ECG. So the
% thoracic is a "noise only" signal that serves as a reference to denoise
% the abdomen signal.

% And the signal/noise ratio differs from one abdomen electrode to another.
% It depends on the position/distance electrode-baby and electrode-mother's
% heart. So if the electrode lies closer to the mother's heart its signal
% will have more noise. It also depends on the environment the wave has
% to cross to get to the electrode.

%% III. Widrow algorithm (LMS)

% Algo inputs:
% - a signal channel: abdomen_ecg
% - one, or more (multi-references version), reference: thorax_ecg
% Algo outputs:
% - maternal ECG
% - ECG of the fetus

% LMS filter settings

L = 8; % RIF filter order associated with the reference signal

bound_mu = 2/(L*mean(thorax_ecg.^2)); % Upper terminal of the adaptivity step in order to guarantee covergence (CS)

% the hypotheses from which the sufficient conditions which ensure the
% cvgce are not necessarily verified ... We can have an approximation an
% then experiment until we have a good convergence.

% In order to have the necessary adaptability and convergence, we take a
% value below the limit for the adaptability step of the LMS algo
% There are 3 cases, with d the denominator: 
mu = bound_mu/1000; mu2 = bound_mu/18; mu3 = bound_mu/10; 

% After some experiments, I found out the choice mu = bound_mu/18 was one
% of the best choices. If we go further to bound_mu/18 the algorithm will diverge. 
% And there's no visible difference between bound_mu/18 and bound_mu/15 anyway.
% The first choice of bound_mu/1000 was also not good because when the
% filter converges that slowly, we do not have enough data for it to finish
% converging. Even between /18 and /100 there is no clear difference.

% In short, if mu is too big, it causes instability in the algorithm
% (fluctuations of coefficients, divergence). If mu is too small, we cannot
% converge with the data we have.

% So we can use this method of experimenting many times combined with the
% mathematical approximations for the best step size mu.

% Declaration / pre-allocation
ecg_maternal_estim = zeros(1, N); ecg_maternal_estim2 = zeros(1, N); ecg_maternal_estim3 = zeros(1, N); % maternal cardio estimation
error = zeros(1, N); error2 = zeros(1, N); error3 = zeros(1, N); % prediction error
w = zeros(L, N); w2 = zeros(L, N); w3 = zeros(L, N); % dimensions 

% Initilization
error(1) = abdomen_ecg(1); error2(1) = abdomen_ecg(1); error3(1) = abdomen_ecg(1);
w(1, 1) = mu*error(1)*thorax_ecg(1); w2(1, 1) = mu2*error2(1)*thorax_ecg(1); w3(1, 1) = mu3*error3(1)*thorax_ecg(1);
ecg_maternal_estim(1) = w(1, 1)*thorax_ecg(1); ecg_maternal_estim2(1) = w2(1, 1)*thorax_ecg(1); ecg_maternal_estim3(1) = w3(1, 1)*thorax_ecg(1);

% Iterations 
for k = 2 : N
    if k < L
        y_k = 0; y2_k = 0; y3_k = 0;
    else
        y_k = thorax_ecg(k - L + 1 : k)'; y2_k = thorax_ecg(k - L + 1 : k)'; y3_k = thorax_ecg(k - L + 1 : k)';
    
        error(k) = abdomen_ecg(k) - w(:, k - 1)'*y_k; error2(k) = abdomen_ecg(k) - w2(:, k - 1)'*y2_k; error3(k) = abdomen_ecg(k) - w3(:, k - 1)'*y3_k; % prediction error 
        w(:, k) = w(:, k - 1) + mu*error(k)*y_k; w2(:, k) = w2(:, k - 1) + mu2*error2(k)*y2_k; w3(:, k) = w3(:, k - 1) + mu3*error3(k)*y3_k; % update filter
        ecg_maternal_estim(k) = w(:, k)'*y_k; ecg_maternal_estim2(k) = w2(:, k)'*y2_k; ecg_maternal_estim3(k) = w3(:, k)'*y3_k; % find maternal ECG at k
    end
end

ecg_fetal_estim = abdomen_ecg - ecg_maternal_estim; ecg_fetal_estim2 = abdomen_ecg - ecg_maternal_estim2; ecg_fetal_estim3 = abdomen_ecg - ecg_maternal_estim3;

%% Visualization of results
% The filter's first coefficients
figure();
subplot(311);
plot(time, w(1, :)'); % change 1 to other numbers to see other coefficients
grid on ; 
xlabel('Time (seconds)');
title('The 1st coeff, d = 1000');
subplot(312);
plot(time, w2(1, :)'); % change 1 to other numbers to see other coefficients
grid on ; 
xlabel('Time (seconds)');
title('The 1st coeff, d = 18');
subplot(313);
plot(time, w3(1, :)'); % change 1 to other numbers to see other coefficients
grid on ; 
xlabel('Time (seconds)');
title('The 1st coeff, d = 10');

% Looking at this we see that the 1st coeff for small mu doesn't have
% enough data to converge. The 2nd one seems the best: it still has some
% small fluctuations but it is able to converge to the right value. The 3rd
% one clearly diverges as the step was too big. For d = 18 we have arrived
% at a minima which gives good convergence so we should stay there.

% So the lesson here is we should always check the coefficients to verify
% if the algorithm works well with the chosen mu or not. We should aim to
% have the fastest convergence possible without much fluctuations so this 
% can be applied in real-time. Clearly the 1st case could not converge 
% within 10 seconds. 

% I verified all the 8 coeff and noticed the same properties for them. A
% lesson is if we see that a coeff keeps increasing or decreasing (like the
% 1st case: each of the 8 coeff either goes up or down all the time), it
% means that the chosen mu is not large enough to converge. The sign of
% divergence is very clear to see. So if we don't have these 2 it means the
% convergence is OK.

%% The error and the signals
figure();
subplot(231);
plot(time, ecg_maternal_estim', 'r', time, thorax_ecg', 'b'); 
grid on ; 
xlabel('Time (seconds)');
ylabel('Signal (voltage)');
legend('Estimated maternal ECG', 'Thorax');
title('The estimated maternal ECG & the thorax, d = 1000');
subplot(232);
plot(time, ecg_maternal_estim2', 'r', time, thorax_ecg', 'b'); 
grid on ; 
xlabel('Time (seconds)');
ylabel('Signal (voltage)');
legend('Estimated maternal ECG', 'Thorax');
title('The estimated maternal ECG & the thorax, d = 18');
subplot(233);
plot(time, ecg_maternal_estim3', 'r', time, thorax_ecg', 'b'); 
grid on ; 
xlabel('Time (seconds)');
ylabel('Signal (voltage)');
legend('Estimated maternal ECG', 'Thorax');
title('The estimated maternal ECG & the thorax, d = 10');
subplot(234);
plot(time, error'); 
grid on ; 
xlabel('Time (seconds)');
ylabel('Signal (voltage)');
title('The error, d = 1000');
subplot(235);
plot(time, error2'); 
grid on ; 
xlabel('Time (seconds)');
ylabel('Signal (voltage)');
title('The error, d = 18');
subplot(236);
plot(time, error3'); 
grid on ; 
xlabel('Time (seconds)');
ylabel('Signal (voltage)');
title('The error, d = 10');

% For the error, it is overall decreased given a large enough mu (like for 
% d = 1000 it went to around 30 but for d = 18 it stayed at around less 
% than 20). 
% It varies around 0 but we have a peak at the estimation of peaks in the 
% maternal ECG. The higher the mu the more frequent the peaks.

%% The signals before and after filtering
figure();
subplot(321);
plot(time, abdomen_ecg'); 
grid on ; 
xlabel('Time (seconds)');
title('The abdomen ECG');
subplot(322);
plot(time, abdomen_ecg'); 
grid on ; 
xlabel('Time (seconds)');
title('The abdomen ECG');
subplot(323);
plot(time, ecg_maternal_estim'); 
grid on ; 
xlabel('Time (seconds)');
title('The estimated maternal ECG, d = 1000');
subplot(324);
plot(time, ecg_maternal_estim2'); 
grid on ; 
xlabel('Time (seconds)');
title('The estimated maternal ECG, d = 18');
subplot(325);
plot(time, ecg_fetal_estim'); 
grid on ; 
xlabel('Time (seconds)');
title('The estimated fetal ECG, d = 1000');
subplot(326);
plot(time, ecg_fetal_estim2'); 
grid on ; 
xlabel('Time (seconds)');
title('The estimated fetal ECG, d = 18');

% This figure helps us see clearer the effect of mu on convergence through 
% the estimation of peaks. For very small mu, the estimation of peaks at 
% the first few seconds was not fine. The quality of estimation increases 
% after some time but the 2nd case it took much less time. The estimation 
% of peaks was correct right at the 1st peak. 

%% Experiment with different order L
% From here on, we will use mu = Bu/18 as we already saw that it give a fast
% convergence.
% I experimented with different values of L by changing it. I put the
% results at the end of this report.

% What I observed: the bigger the L, the smaller the step size mu allowed,
% but if the mu is small enough the convergence should still be guaranteed. 
% So the algorithm has more 0 at the beginning and we lose some
% information. But we will have more weights so more filtering (the filter
% is more complex). On the other hand if we don't have enough weights we 
% won't be able to capture information with high variations (in DSP). The
% algorithm will diverge if L is too low.

% However we need to store and perform calculations on the weights to this 
% would consume memory. There's no point in increasing the complexity if
% the results are already good enough in terms of visibility. We stop at
% the point where further increase in complexity gives almost no further 
% increase in the clarity of results.

% In this case L = 8 is OK.

%% IV. Multi-reference Widrow algorithm (LMS)
% In this interesting part we will filter the signal using 3 "nois only" 
% references. All the 3 thoracic signals are taken, giving us a 
% multi-reference filter.

load foetal_ecg.mat;   % load data
Fs = 250 ; % sampling frequency in Hz
N = size(foetal_ecg, 2) ; % number of samples
time = (0 : N - 1)/Fs ; % discrete time

ireferences = [6 7 8]; % the reference electrode (s) / signals
nRef = length(ireferences) ; % number of reference signals
thorax_ecg_mult = foetal_ecg(ireferences, 1 : N) ; % nRef x N <- multi-ref

% Signal to denoise
elec_abdo = 1;
abdomen_ecg = foetal_ecg(elec_abdo, :); 

clear foetal_ecg ;

L_mult = 8;

% choose the smallest bound which corresponds to 
bound_mu_mult = 2/(L_mult*max(mean(thorax_ecg_mult.^2))); 
% Overall the multi-ref gives us a wider choice of mu.

mu_mult = bound_mu_mult/2;
% Overall the multi-ref gives us a wider choice of mu.

ecg_maternal_estim_mult = zeros(1, N); % maternal cardio estimation
error_mult = zeros(1, N); % prediction error
w_mult = zeros(L_mult, length(ireferences), N); % dimensions 

% There is now a response by reference electrode
% The LMS filter being the sum of each of these "nRef" filters
% The responses of the (sub) filters at time k will be stored in the 
% vectors w_mult (:, 1, k), w_mult (:, 2, k), ...

% Initilization
error_mult(1) = abdomen_ecg(1);
w_mult(1, :, 1) = mu*error_mult(1)*thorax_ecg_mult(:, 1)';
ecg_maternal_estim_mult(1) = w_mult(1, :, 1)*thorax_ecg_mult(:, 1); 
     
% Prediction error: the prediction is now the sum of the predictions obtained 
% for each of the nRef = 3 abdominal electrodes using subdose responses are 
% stored in w_mult(:, 1, k-1), w_mult(:, 2, k-1) and w_mult(:, 3, k-1)
    
for k = 2 : N
    if k < L
        y_k = 0; 
    else
        y_k = thorax_ecg_mult(:, k - L + 1 : k)';
    
        sum = 0;
        for j = 1 : length(ireferences)
            sum = sum + w_mult(:, j, k - 1)'*y_k(:, j);
        end    
    
        error_mult(k) = abdomen_ecg(k) - sum; % prediction error 
        w_mult(:, :, k) =  w_mult(:, :, k - 1) + mu_mult*error_mult(k)*y_k; % update filter
    
        ecg_maternal_estim_mult(k) = 0;
        for i = 1 : length(ireferences)
            ecg_maternal_estim_mult(k) = ecg_maternal_estim_mult(k) + w_mult(:, i, k)'*y_k(:, i);
        end
    end
end

ecg_fetal_estim_mult = abdomen_ecg - ecg_maternal_estim_mult;

% Comparison with the single-reference case

figure();
subplot(321);
plot(time, abdomen_ecg'); 
grid on ; 
xlabel('Time (seconds)');
xlim([2 4]);
title('The abdomen ECG');
subplot(322);
plot(time, abdomen_ecg'); 
grid on ; 
xlabel('Time (seconds)');
xlim([2 4]);
title('The abdomen ECG');
subplot(323);
plot(time, ecg_maternal_estim2'); 
grid on ; 
xlabel('Time (seconds)');
xlim([2 4]);
title('The estimated maternal ECG, single-reference');
subplot(324);
plot(time, ecg_maternal_estim_mult'); 
grid on ; 
xlabel('Time (seconds)');
xlim([2 4]);
title('The estimated maternal ECG, multi-reference');
subplot(325);
plot(time, ecg_fetal_estim2'); 
grid on ; 
xlabel('Time (seconds)');
xlim([2 4]);
title('The estimated fetal ECG, single-reference');
subplot(326);
plot(time, ecg_fetal_estim_mult'); 
grid on ; 
xlabel('Time (seconds)');
xlim([2 4]);
title('The estimated fetal ECG, multi-reference');

% If we look closely, especially the 2nd figure each column, we see that the
% estimation is much better with multi-ref: the estimated maternal ECG is
% smoother which means it contains only the mother's ECG. If we look in the
% results, where the maternal ECG used to be, we'll see that the multi-ref
% case has much fewer peaks (or variations). So the signal was better
% filtered and we can better detect the heartbeat of the foetus. 
% This can be explained by the fact that here the multiple 
% references can 'help' each other through the cross-effect among them.

% After all, the objective of this is to help the doctors measure the
% heartbeat of the foetus, so having multiple thoracic electrodes will enable
% the multi-ref algorithm and improve visual the quality.


%% V. Conclusion
% After this Lab, I understood the use of the Widrow filter to denoise
% signals, especially in the case the noise is bigger in amplitude than the
% useful signal and we can measure the "noise only" reference. Besides, the
% version with multi-reference gaves smoother results than the one with a
% single reference.

% Below I put the extra figures for demonstration.
