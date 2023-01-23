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
Pb = struct('T1',770e-3,'T2',33e-3,'f',1e-3,'k',200,'dw',700);
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
clc;clear; close all
Pa = struct('T1',1,'T2',0.2); % s
Pb = struct('T1',0.1,'T2',0.1,'f',0.0003636,'k',500,'dw',20000);
sys = struct('w1',512,'tp',10,'offsets',30e3:-100:-10e3);
Z = CEST_multipool(sys,Pa,Pb);
plot(sys.offsets,Z); set(gca,'XDir','reverse')
%% Alter kb
clc;clear;
kb = [300 500 900 3000 7000 9000]; % Hz
Pa = struct('T1',0.5,'T2',0.2); % s
Pb = struct('T1',0.1,'T2',0.1,'f',0.0003636,'dw',20000);
sys = struct('w1',512,'tp',10,'offsets',30e3:-100:-10e3);
colors = lines(length(kb));
h = zeros(length(kb),1);
for i=1:length(kb)
    Pb.k = kb(i);
    plot(sys.offsets,CEST_multipool(sys,Pa,Pb),'color',colors(i,:),'DisplayName',strcat('kb= ',num2str(kb(i)),' Hz'));
    hold on
end
set(gca,'XDir','reverse')
legend