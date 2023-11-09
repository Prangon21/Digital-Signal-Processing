clc;
clear all;
close all;

%% Message Signal Generation 

[message1, fs] = audioread('message_1.wav');
[message2, fs] = audioread('message_2.wav');

m1 = message1(:,1)'  
m2 = message2(:,1)'  

ml1=length(m1);
ml2=length(m2);

n  = ceil(-(ml1)/2):floor((ml1-1)/2); % since ml1 = ml2
ts = 1/fs;
t  = n*ts;

N  = 2.^nextpow2(ml1);
fn = [0:1/N:1-1/N]*fs-fs/2; % Frequency axis for spectrum
tn = 1./fn;

figure(1)
subplot(211)
plot(t,message1);
title("Audio Input 1; m1(t)");xlabel("Time");ylabel("Magnitude");

subplot(212)
plot(t,message2);
title("Audio Input 2; m2(t)");xlabel("Time");ylabel("Magnitude");

%% Frequency Spectrum of Message Signal

mw1 = fft(m1,N);
MW1 = mw1/fs; 

mw2 = fft(m2,N);
MW2 = mw2/fs; 

figure(2)
subplot(211)
plot(fn, abs(fftshift(MW1)));
subplot(212)
plot(fn, abs(fftshift(MW2)));

%% Low pass or Anti-aliasing Filter Output

fcut = 4000;    % Cut-off frequency of the low pass filter
ncut = 4000;
H = zeros(1,N); % Filter vector zeros(row,column)
H(1:ncut) = ones(1,ncut) % Low pass filter with gain 1
H(N-ncut+1:N) = ones(1, ncut); % Other portion of the low pass filter

Ufiltered1 = real(MW1).*H; % Filtering the modulated spectrum
Ufiltered2 = real(MW2).*H;

figure(3)

subplot(211)
plot(fn, abs(fftshift(Ufiltered1)));
title('Spectrum of the filtered Signal');
subplot(212)
plot(fn, abs(fftshift(Ufiltered2)));
title('Spectrum of the filtered Signal');

%% IFFT OF FILTERED SIGNAL

ufiltered_t_1 = ifft(Ufiltered1)*fs;
ufiltered_t_2 = ifft(Ufiltered2)*fs;

%% Carrier Signal Generation

fc1 = 1140630120e3;
fc2 = 1140630120;  % Dividing with 3 to remove the overlapping
c1 = cos(2*pi*fc1*t); % Carrier signal 1
c2 = cos(2*pi*fc2*t); % Carrier signal 2

figure(10)

subplot(211);plot(t,c1);
subplot(212);plot(t,c2);

C1 = fft(c1,N);
C1 = C1/fs;
carrier_1 = real(ifft(C1))*fs; % for array matching

C2 = fft(c2,N);
C2 = C1/fs; 
carrier_2 = real(ifft(C2))*fs; % for array matching

%% Modulation

modulated_1 = ufiltered_t_1.*carrier_1;
modulated_2 = ufiltered_t_2.*carrier_2;

MODULATED_1 = fft(modulated_1,N);
MODULATED_1 = MODULATED_1/fs;

MODULATED_2 = fft(modulated_2,N);
MODULATED_2 = MODULATED_2/fs;

figure(4)
subplot(211)
plot(fn,fftshift(abs(MODULATED_1)))
title(' Modulated Audio Data 1 ')

subplot(212)
plot(fn,fftshift(abs(MODULATED_2)))
title(' Modulated Audio Data 2 ')

%% Adder Operation

Added_Modulated = real(MODULATED_1) + real(MODULATED_2);

ADDED_MODULATED = fft(Added_Modulated,N);
ADDED_MODULATED = ADDED_MODULATED/fs;


figure(5)
plot(fn,fftshift(abs(ADDED_MODULATED)))
title('Added the Modulated Signal in Frequency Domain')

%% Bandpass Filtering

fcut1 = 6000;    % Cut-off frequency of the low pass filter
ncut1 = 4000;
H = zeros(1,N); % Filter vector zeros(row,column)
H(1:ncut) = ones(1,ncut) % Low pass filter with gain 1
H(N-ncut+1:N) = ones(1, ncut); % Other portion of the low pass filter

ADDMODfiltered1 = ADDED_MODULATED.*H; % Filtering the modulated spectrum

figure(6)

subplot(211)
plot(fn, abs(fftshift(ADDMODfiltered1)));
title('Spectrum of the filtered Signal');

fcut1 = 2500;    % Cut-off frequency of the low pass filter
ncut1 = 1500;
H = zeros(1,N); % Filter vector zeros(row,column)
H(1:ncut) = ones(1,ncut) % Low pass filter with gain 1
H(N-ncut+1:N) = ones(1, ncut); % Other portion of the low pass filter

ADDMODfiltered2 = ADDED_MODULATED.*H; % Filtering the modulated spectrum

subplot(212)
plot(fn, abs(fftshift(ADDMODfiltered2)));
title('Spectrum of the filtered Signal');

%% Again ifft

bp_filter_output1 = real(ifft(ADDMODfiltered1))*fs;
bp_filter_output2 = real(ifft(ADDMODfiltered2))*fs;

%% Demodulation

% demod1 = ADDMODfiltered1(:,1)'  
% demod2 = ADDMODfiltered2(:,1)'  

demodulated_1 = bp_filter_output1.*carrier_1;
demodulated_2 = bp_filter_output2.*carrier_2;

DEMODULATED_1 = fft(demodulated_1,N);
DEMODULATED_1 = DEMODULATED_1/fs;

DEMODULATED_2 = fft(demodulated_2,N);
DEMODULATED_2 = DEMODULATED_2/fs;

figure(7)
subplot(211)
plot(fn,fftshift(abs(DEMODULATED_1)))
title(' Demodulated Audio Data 1 ')

subplot(212)
plot(fn,fftshift(abs(DEMODULATED_2)))
title(' Demodulated Audio Data 2 ')

%% Low pass or Reconstruction Filter Output

Lastfiltered1 = DEMODULATED_1.*H; % Filtering the modulated spectrum
Lastfiltered2 = DEMODULATED_2.*H;

figure(8)

subplot(211)
plot(fn, abs(fftshift(Lastfiltered1)));
title('Spectrum of the filtered Signal');
subplot(212)
plot(fn, abs(fftshift(Lastfiltered2)));
title('Spectrum of the filtered Signal');

figure(9)

subplot(211)
plot(t, real(ifft(M1))*fs);
title('Spectrum of the filtered Signal');
subplot(212)
plot(t, real(ifft(M2))*fs);
title('Spectrum of the filtered Signal');

