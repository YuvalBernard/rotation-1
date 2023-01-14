clc;clear;close all;

Pa = struct('T1',3,'T2',2);
Pb = struct('T1',770e-3,'T2',33e-3,'f',1e-3,'k',200,'dw',-700);
sys = struct('w1',3*42.58,'offsets',1000:-10:-1000,'tp',10);

[Z,MTR] = CEST_multipool(sys,Pa,Pb);
plot(sys.offsets, Z,'r-'); ax = gca; ax.XDir = 'reverse';
hold on
Pb.f = 0;
Z_ref = CEST_multipool(sys,Pa,Pb);
plot(sys.offsets, Z_ref,'b-');
plot(sys.offsets,Z_ref-Z,'m--')
plot(sys.offsets,MTR,'c.')

%% Change kb
clc;clear;close all;

Pa = struct('T1',3,'T2',2);
Pb = struct('T1',770e-3,'T2',33e-3,'f',0.2,'dw',-2500);
sys = struct('w1',500,'offsets',5000:-10:-5000);
k = [10 50 200 500 1000];
hold on; h = zeros(size(k));
for i = 1:length(k)
    Pb.k = k(i);
    Z = CEST_multipool(sys,Pa,Pb);
    h(i) = plot(sys.offsets,Z,'DisplayName',['k= ',num2str(k(i)),' Hz']);
end
ax = gca; ax.XDir = 'reverse'; legend(h(1:end))

%% Change w1
clc;clear;close all;

Pa = struct('T1',3,'T2',2);
Pb = struct('T1',770e-3,'T2',33e-3,'f',1e-3,'k',200,'dw',-700);
sys = struct('offsets',1000:-10:-1000,'tp',10);
w1 = [1 2 3]*42.58; % Hz
h = zeros(length(w1),1);
hold on
for i=1:length(w1)
    sys.w1 = w1(i);
    [Z,MTR_asymm] = CEST_multipool(sys,Pa,Pb);
    yyaxis left
    ylabel('Z Specta')
    ylim([-0.2 1])
    h(i) = plot(sys.offsets,Z,'DisplayName',strcat('w1= ',num2str(w1(i)),' Hz'));
    yyaxis right
    ylabel('MTR_{Asymm} Spectra')
    ylim([-0.1 0.5])
    plot(sys.offsets,MTR_asymm);
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