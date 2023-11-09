clc;
clear all;
close all;

% DSB-SC Modulation and demodulation of message signal using FFT

%% message signal and carrier signal generation

fc = 1000; % Carrier frequency
[Y, fs]=audioread('message_1.wav');

m=Y(:,1)';% 1 ch message signal 
% Make a carrier vector of length equal to the input vector m 
ml=length(m);
n=ceil(-(ml)/2):floor((ml-1)/2);
ts=1/fs;
t=n*ts;
c = cos(2*pi*fc*t); % carrier signal

%% Modulation
u = m.*c; % DSB-AM modulated signal
% N = 2^15; % FFT Bin size
N = 2.^nextpow2(length(u)); % FFT Bin size

M = fft(m, N); % N point DFT of message signal
% M = M/fs; % scaling
C = fft(c, N); % N point DFT of message signal
% M = M/fs; % scaling

U = fft(u, N); % N point DFT of modulated signal
% U = U/fs;

% Visualizing modulated spectrum
fn = [0:1/N:1-1/N]*fs-fs/2; % Frequency axis for spectrum
figure(1)
subplot(311), plot(fn, abs(fftshift(M)), 'Linewidth', 2); % showing message spectrum
title('Spectrum of the message signal');
subplot(312), plot(fn, abs(fftshift(C)), 'Linewidth', 2); % showing message spectrum
title('Spectrum of the message signal');
subplot(313), plot(fn, abs(fftshift(U)), 'Linewidth', 2); % Showing the modulated spectrum
title('Spectrum of the modulated signal');


