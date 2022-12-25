clc;clear;close all;

T1a = 3; T2a = 2; % s
T1b = 770e-3; T2b = 33e-3; % s
M0a = 1; M0b = M0a*1e-3; % arb
kb = 200; % s^-1
w1 = 42.58*3; % Hz (irradiation intensity)
dwa = 1000:-10:-1000; % Hz (Saturation at given freqs)
db = -700; % offset between wa and wb; dwb = dwa + db
[Z,A,domain] = CEST(T1a,T2a,T1b,T2b,kb,M0a,M0b,dwa,db,w1,10);
plot(dwa,Z,'k',dwa(1:domain),A,'b'); ax = gca; ax.XDir = 'reverse';
%% Change w1
clc;clear;close all;

T1a = 3; T2a = 2; % s
T1b = 770e-3; T2b = 33e-3; % s
M0a = 1; M0b = M0a*1e-3; % arb
kb = 200; % Hz
dwa = 1000:-10:-1000; % Hz (Saturation at given freqs)
db = -700; % offset between wa and wb; dwb = dwa + db
w1 = [1 2 3]*42.58; % Hz
h = zeros(length(w1),1);
hold on
for i=1:length(w1)
    [Z,A,domain] = CEST_krylov(T1a,T2a,T1b,T2b,kb,M0a,M0b,dwa,db,w1(i),10);
    yyaxis left
    ylabel('Z Specta')
    ylim([-0.2 1])
    h(i) = plot(dwa,Z,'DisplayName',strcat('w1= ',num2str(w1(i)),' Hz'));
    yyaxis right
    ylabel('A Spectra')
    ylim([-0.1 0.5])
    plot(dwa(1:domain),A);
end
ax = gca;
ax.XDir = 'reverse';
legend(h(1:end))
%% Replicate experimental results
% Fit data from "Numerical Solution of the Bloch Equations Provides
% Insights Into the Optimum Design of PARACEST Agents
% for MRI".
% Goal: Obtain similar results to Fig 1.
clc;clear;close all;
kb = 500; % Hz
dwa = 30000:-100:-10000; % Hz
db = -20000;
w1 = 512; % Hz
T1a = 1; T2a = 0.2; % s
T1b = 0.1; T2b = 0.1; % s
M0a = 1; M0b = 0.0003636;
[Z,A,domain] = CEST(T1a,T2a,T1b,T2b,kb,M0a,M0b,dwa,db,w1,10);

plot(dwa,Z,dwa(1:domain),A)
ax = gca;
ax.XDir = 'reverse';
%% Alter kb
clc;clear;close all;
kb = [300 500 900 3000 7000 9000]; % Hz
dwa = 30000:-100:-10000; % Hz
db = -20000;
w1 = 512; % Hz
T1a = 2; T2a = 0.2; % s
T1b = 0.1; T2b = 0.1; % s
M0a = 1; M0b = 0.0003636;
colors = lines(length(kb));
h = zeros(length(kb),1);
for i=1:length(kb)
    Z = CEST(T1a,T2a,T1b,T2b,kb(i),M0a,M0b,dwa,db,w1);
    h(i) = plot(dwa,Z,'color',colors(i,:),'DisplayName',strcat('kb= ',num2str(kb(i)),' Hz'));
    hold on
end
ax = gca;
ax.XDir = 'reverse';
legend(h(1:end))