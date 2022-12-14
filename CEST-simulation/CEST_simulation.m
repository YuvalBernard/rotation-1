clc;clear;close all;

T1a = 3; T2a = 2; % s
T1b = 770e-3; T2b = 33e-3; % s
M0a = 1; M0b = M0a*5e-3; % arb
kb = 200; % s^-1
w1 = 127.7; % Hz (irradiation intensity)
dwa = -1000:10:1000; % Hz (Saturation at given freqs
dwb = dwa + 700;

[Z,A,domain] = CEST(T1a,T2a,T1b,T2b,kb,M0a,M0b,dwa,dwb,w1);

plot(dwa,Z,'k--',dwa(1:domain),A,'b')

%% Change w1
w1 = [42.6 85.2 127.7];
colors = lines(length(w1));
h = zeros(length(w1),1);
for i=1:length(w1)
    [Z,A,domain] = CEST(T1a,T2a,T1b,T2b,kb,M0a,M0b,dwa,dwb,w1(i));
    h(i) = plot(dwa,Z,'color',colors(i,:),'LineStyle','--','DisplayName',strcat('w1= ',num2str(w1(i)),' Hz'));
    hold on
    plot(dwa(1:domain),A,'color',colors(i,:),'LineStyle','-');
end
legend(h(1:end))
%% Replicate experimental results
% Fit data from "Numerical Solution of the Bloch Equations Provides
% Insights Into the Optimum Design of PARACEST Agents
% for MRI".
% Goal: Obtain similar results to Fig 1.
clc;clear;close all;
kb = 3000; % Hz
dwa = -10000:100:30000; % Hz
dwb = dwa + 20000;
w1 = 2500; % Hz
T1a = 2; T2a = 0.2; % s
T1b = 0.1; T2b = 0.08; % s
M0a = 1; M0b = 0.0007272*M0a;
[Z,A,domain] = CEST(T1a,T2a,T1b,T2b,kb,M0a,M0b,dwa,dwb,w1);

plot(dwa,Z)