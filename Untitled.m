
clc;
clear all;
close all;


Fs = 16000;       
dt = 1/Fs;        

F = 10e3;                  
D = .8;                       
PW = D*F;                  
f = 1/F;                   
t = -T/2:dt:T/2;             
n = t/dt;                    
L = PW/dt;                    
x = zeros(1,length(t));  



x(find(abs(n)<=L/2))=1;  
plot(t,x, 'Linewidth', 2);