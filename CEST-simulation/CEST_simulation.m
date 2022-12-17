clc;clear;close all;

T1a = 3; T2a = 2; % s
T1b = 770e-3; T2b = 33e-3; % s
M0a = 1; M0b = M0a*1e-3; % arb
kb = 200; % s^-1
w1 = 127.7; % Hz (irradiation intensity)
dwa = 1000:-0.1:-1000; % Hz (Saturation at given freqs)
db = -700; % offset between wa and wb; dwb = dwa + db
[Z,A,domain] = CEST(T1a,T2a,T1b,T2b,kb,M0a,M0b,dwa,db,w1);
plot(dwa,Z,'k',dwa(1:domain),A,'b');
ax = gca;
ax.XDir = 'reverse';
%% Change w1
clc;clear;close all;

T1a = 3; T2a = 2; % s
T1b = 770e-3; T2b = 33e-3; % s
M0a = 1; M0b = M0a*1e-3; % arb
kb = 200; % s^-1
dwa = 1000:-0.1:-1000; % Hz (Saturation at given freqs)
db = -700; % offset between wa and wb; dwb = dwa + db
w1 = [42.6 85.2 127.7];
colors = lines(length(w1));
h = zeros(length(w1),1);
for i=1:length(w1)
    [Z,A,domain] = CEST(T1a,T2a,T1b,T2b,kb,M0a,M0b,dwa,db,w1(i));
    h(i) = plot(dwa,Z,'color',colors(i,:),'DisplayName',strcat('w1= ',num2str(w1(i)),' Hz'));
    hold on
    plot(dwa(1:domain),A,'color',colors(i,:),'LineStyle','-');
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
[Z,A,domain] = CEST(T1a,T2a,T1b,T2b,kb,M0a,M0b,dwa,db,w1);

plot(dwa,Z,dwa(1:domain),A)
ax = gca;
ax.XDir = 'reverse';