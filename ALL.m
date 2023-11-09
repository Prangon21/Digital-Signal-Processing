clc;
clear all;
close all;


%% Audio Inputs in Time Domain

fs = 16000;

m1 = audioread('message_1.wav');
m2 = audioread('message_2.wav');

message1 = m1.';
message2 = m2.';

t1 = linspace(0,length(message1)/fs,length(message1));
t2 = linspace(0,length(message2)/fs,length(message2));

figure(1)
subplot(211)
plot(t1,message1);
title("Audio Input 1; m1(t)");xlabel("Time");ylabel("Magnitude");grid on

subplot(212)
plot(t2,message2);
title("Audio Input 2; m2(t)");xlabel("Time");ylabel("Magnitude");grid on

%% Audio Inputs in Frequency Domain

N = 2.^nextpow2(length(message1));
fn = [0:1/N:1-1/N]*fs-fs/2; % Frequency axis for spectrum

mm1 = fft(m1,N);plot(abs(mm1));
MM1 = mm1/fs; 
mm2 = fft(m2,N);
MM2 = mm2/fs; 

figure(2)
subplot(211)
plot(fn, abs(fftshift(MM1)))
subplot(212)
plot(fn, abs(fftshift(MM2)))

%% Filtering

fcut = 4e3; % cut-off frequency of the low pass filter
ncut = floor(fcut*fs/N); % index vector for filter upto cut-off frequency

H = zeros(1,N); % Filter vector
H(1:ncut) = 1*ones(1,ncut);      % Low pass filter with gain 1
H(N-ncut+1:N) = 1*ones(1, ncut); % Other portion of the low pass filter
Ufiltered1 = MM1.*H; % Filtering the modulated spectrum
Ufiltered2 = MM2.*H;

figure(3)
subplot(311)
plot(H, 'linewidth',2);
subplot(312)
plot(fn, abs(fftshift(Ufiltered1)), 'Linewidth', 2);
title('Spectrum of the filtered Signal');
subplot(313)
plot(fn, abs(fftshift(Ufiltered2)), 'Linewidth', 2);
title('Spectrum of the filtered Signal');
