% CEST and DEST simulations (Slides from 19.12.22)
clc;clear;close all;

%% Lithium CEST
% Interaction between Li dendrites (free) and SEI (bound)
% Experiment done on Li-7 spins at 9.4T
gamma = 16.546; % MHz/T; gyromagnetic ratio
B0 = 9.4; % T

T1a = 100; T2a = 1/(pi*40); % s
T1b = 150e-3; T2b = 0.5e-3; % s
M0a = 1; M0b = M0a*5e-3; % arb
kb = linspace(10,1000,5); % s^-1 Select values between 10-1000 Hz
w1 = 2000; % Hz (irradiation intensity)
dwa = 500*B0*gamma:-1000:-500*B0*gamma; % Hz (Saturation at given freqs)
db = -260*B0*gamma; % offset between wa and wb; dwb = dwa + db
h = zeros(length(w1),1);
figure
for i=1:length(kb)
    [Z,A,domain] = CEST(T1a,T2a,T1b,T2b,kb(i),M0a,M0b,dwa,db,w1);
    subplot(2,1,1)
    hold on
    ylabel('Z Spectra');  ax = gca; ax.XDir = 'reverse';
    h(i) = plot(dwa,Z,'DisplayName',strcat('kb= ',num2str(kb(i)),' Hz'));
    legend(h(1:end))
    subplot(2,1,2)
    hold on; xlim([dwa(end) dwa(1)])
    ylabel('A Spectra'); ax = gca; ax.XDir = 'reverse';
    plot(dwa(1:domain),A)
end