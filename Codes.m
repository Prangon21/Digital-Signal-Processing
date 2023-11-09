clear all;
close all;
clc;

%% Audio Inputs in Time Domain

Fs = 16e3;

m1 = audioread('message_1.wav');
m2 = audioread('message_2.wav');
t1 = linspace(0,length(m1)/Fs,length(m1));
t2 = linspace(0,length(m2)/Fs,length(m2));

figure(1)
subplot(211)
plot(t1,m1);
title("Audio Input 1; m1(t)");xlabel("Time");ylabel("Magnitude");grid on
xlim([0,11]);

subplot(212)
plot(t2,m2);
title("Audio Input 2; m2(t)");xlabel("Time");ylabel("Magnitude");grid on
xlim([0,11]);

%% Audio Inputs in Frequency Domain

N  = pow(3 ; % FFT point number
fn = [0:1/N:1-1/N]*Fs-Fs/2; % Frequency axis for spectrum

mm1 = fft(m1,N);plot(abs(mm1));
MM1 = mm1/Fs; 
mm2 = fft(m2,N);
MM2 = mm2/Fs; 

figure(2)
subplot(211)
plot(fn, abs(fftshift(MM1)))
subplot(212)
plot(fn, abs(fftshift(MM2)))

%% Filtering

% fcut = 4e3; % cut-off frequency of the low pass filter
% ncut = floor(fcut*Fs/N); % index vector for filter upto cut-off frequency
% H = zeros(1,N); % Filter vector
% H(1:ncut) = 1*ones(1,ncut);      % Low pass filter with gain 1
% H(N-ncut+1:N) = 1*ones(1, ncut); % Other portion of the low pass filter
% Ufiltered1 = MM1.*H; % Filtering the modulated spectrum
% Ufiltered2 = MM2.*H;
% 
% figure(3)
% subplot(311)
% plot(H, 'linewidth',2);
% subplot(312)
% plot(fn, abs(fftshift(Ufiltered1)), 'Linewidth', 2);
% title('Spectrum of the filtered Signal');
% subplot(313)
% plot(fn, abs(fftshift(Ufiltered2)), 'Linewidth', 2);
% title('Spectrum of the filtered Signal');

%% Modulation

% fc1 = 1140630120e3
% fc2 = fc1*6;
% carrier_1 = cos(2*pi*fc1*t); % carrier signal 1
% carrier_2 = cos(2*pi*fc2*t); % carrier signal 2
% 
% modulated_1 = m1.*carrier_1;
% modulated_2 = m2.*carrier_2;
% 
% MODULATED_1 = fft(modulated_1,N);
% MODULATED_1 = MODULATED_1/Fs;
% 
% MODULATED_2 = fft(modulated_2,N);
% MODULATED_2 = MODULATED_2/Fs;
% 
% Added_Modulated = modulated_1 + modulated_2;
% ADDED_MODULATED = fft(Added_Modulated,N)
% ADDED_MODULATED = Added_Modulated/Fs;


fc_1 = 1140630120e3 ;   % Carrier Frequency 1 = Sum of all the digits of your group members = 180105028 + 180105035 + 180105040 + 180105042 + 180105044 = 900525189
fc_2 = fc_1*6 ;         % Carrier Frequency 2

Len = length(m1) ;
F_axis = fc_1*20 ;      % Defining F axis
fs = F_axis ; 
ts = 1/fs ;
F = (0 : 1/Len : 1-(1/Len))*fs - (fs/2) ;

A = length(m1)/2 ;
ty = -A*ts : ts : A*ts-ts ;     % defining y axis

Carrier_1 = cos(2*pi*fc_1*ty) ; % Carrier signal 1
Carrier_2 = cos(2*pi*fc_2*ty) ; % Carrier signal 2

C_1 = fft(Carrier_1) ;
C_2 = fft(Carrier_2) ;

figure(4)
subplot(323)
plot(ty, Carrier_1)
title(' carrier 1 - time domain ')
grid on
xlim([0,4e-9]);

subplot(324)
plot(F,fftshift(abs(C_1)))
title(' carrier 1 - frequency domain ')
grid on
xlim([-1e13,1e13]);

subplot(325)
plot(ty,Carrier_2)
title(' carrier 2 time domain ')
grid on
xlim([0,4e-9]);

subplot(326)
plot(F,fftshift(abs(C_2)))
title(' carrier 2 - frequency domain ')
grid on
xlim([-1e13,1e13]);

subplot(311)
plot(F,fftshift(abs(C_1)))
hold on
grid on
plot(F,fftshift(abs(C_2)))
xlim([-1e13,1e13]);
title(' carrier - frequency domain')






















