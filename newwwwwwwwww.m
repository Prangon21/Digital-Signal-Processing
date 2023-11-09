clc;
clear all;
close all;


%% Message signal and Carrier signal

fc = 1000; % Carrier frequency
[Y, fs]=audioread('message_1.wav');

m=Y(:,1)'  % 1 channel message signal 
ml=length(m);
n=ceil(-(ml)/2):floor((ml-1)/2);
ts=1/fs;
t=n*ts;
c = cos(2*pi*fc*t); % carrier signal

N = 2.^nextpow2(length(m));
fn = [0:1/N:1-1/N]*fs-fs/2; % Frequency axis for spectrum

mm1 = fft(m,N);
MM1 = mm1/fs; 
% mm2 = fft(m2,N);
% MM2 = mm2/fs; 
real = abs(fftshift(MM1));
figure(2)
plot(fn,real)
% subplot(212)
% plot(fn, abs(fftshift(MM2)))

%% Filtering

fcut = 4000; % cut-off frequency of the low pass filter

ncut = 4000 % index vector for filter upto cut-off frequency
H = zeros(1,N); % Filter vector zeros(row,column)
H(1:ncut) = ones(1,ncut) % Low pass filter with gain 1
H(N-ncut+1:N) = ones(1, ncut); % Other portion of the low pass filter

Ufiltered = MM1.*H; % Filtering the modulated spectrum

figure(3)
subplot(211)
plot(H);
xlim([-8000,8000])

subplot(212)
plot(fn, abs(fftshift(Ufiltered)));
title('Spectrum of the filtered Signal');


figure(4)
plot(abs(MM1))

figure(5)
plot(H); 