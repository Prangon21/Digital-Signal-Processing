clc;
clear all;
close all;

Fs = 16000;           % Sampling frequency in hertz
ch = 1;               % Number of channels--2 options--1 (mono) or 2 (stereo)
nbits = 16;           % 8,16,or 24
Nseconds = 10;

for A = 1:1:2;
    if A == 1
        disp('Message 1...')
        disp('Start speaking..')
        recorder1 = audiorecorder(Fs,nbits,ch);
        recordblocking(recorder1,Nseconds);
        disp('End of Recording.');
        x = getaudiodata(recorder1);
        audiowrite('message_1.wav',x,Fs);

    elseif A == 2
        disp('Message 2...')
        recorder2 = audiorecorder(Fs,nbits,ch);
        disp('Start speaking..')
        recordblocking(recorder2,Nseconds);
        disp('End of Recording.');
        y  = getaudiodata(recorder2);
        audiowrite('message_2.wav',y,Fs);
    else disp('Input 1 or 2..')
        
    end
end
