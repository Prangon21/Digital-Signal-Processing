clear all;
close all;
clc;

%% Audio Inputs in Time Domain

Fs = 16e3;

m1 = audioread('message_1.wav');
m2 = audioread('message_2.wav');

% m_1=m1(:,1)';% 1 ch message signal 
% m_2=m2(:,1)';
% Make a carrier vector of length equal to the input vector m 
% 
% message1=length(m_1);
% message2=length(m_2);
% 
% n1=ceil(-(message1)/2):floor((message1-1)/2);
% n2=ceil(-(message2)/2):floor((message2-1)/2);
% 
% ts = 1/fs;
% t1  = n1*ts;
% t2  = n2*ts;

figure(1)
subplot(211)
plot(t1,message1);
title("Audio Input 1; m1(t)");xlabel("Time");ylabel("Magnitude");grid on
% xlim([0,11]);

subplot(212)
plot(t2,message2);
title("Audio Input 2; m2(t)");xlabel("Time");ylabel("Magnitude");grid on
% xlim([0,11]);