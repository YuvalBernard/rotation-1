clc;clear;close all;

T1a = 3; T2a = 2; % s
T1b = 770e-3; T2b = 33e-3; % s
M0a = 1; M0b = M0a*5e-3; % arb
kb = 200; % s^-1
w1 = 127.7; % Hz (irradiation intensity)
dwa = -1000:10:1000; % Hz (Saturation at given freqs
dwb = dwa + 700;
t0 = 0; tmax = 30; q = 1000;

[Z,A,domain] = CEST(T1a,T2a,T1b,T2b,kb,M0a,M0b,dwa,dwb,w1,t0,tmax,q);

plot(dwa,Z,'k--',dwa(1:domain),A,'b')

%% Change w1
for w1 = [42.6 85.2 127.7]
    [Z,A,domain] = CEST(T1a,T2a,T1b,T2b,kb,M0a,M0b,dwa,dwb,w1,t0,tmax,q);
    plot(dwa,Z,'--',dwa(1:domain),A,'-'); hold on
end