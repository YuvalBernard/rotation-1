% CEST and DEST simulations (Slides from 19.12.22)

%% Lithium CEST
clc;clear;close all;
% Interaction between Li dendrites (free) and SEI (bound)
% Experiment done on Li-7 spins at 9.4T
gamma = 16.546; % MHz/T; gyromagnetic ratio
B0 = 9.4; % T

% System parameters:
T1a = [10; 100]; T2a = 1/(pi*35); % s
T1b = 150e-3; T2b = 0.5e-3; % s
M0a = 1; M0b = M0a*5e-3; % arb
kb = [10;50;100;500;1000]; % s^-1 Select values between 10-1000 Hz
w1 = 750; % Hz (irradiation intensity)
dwa = 2000*B0*gamma:-B0*gamma:-2000*B0*gamma; % Hz (Saturation at given freqs)
db = -260*B0*gamma; % offset between wa and wb; dwb = dwa + db

% Calculate spectra and plot results
h = zeros(length(w1),1);
figure
sgtitle('Lithium CEST: Dendrites + SEI')
for j=1:length(T1a)
    for i=1:length(kb)
        [Z,A,domain] = CEST(T1a(j),T2a,T1b,T2b,kb(i),M0a,M0b,dwa,db,w1);
        subplot(2,2,j)
        xline([500 -500]*B0*gamma,'--'); xlabel('ppm')
        title(strcat('T1a=',num2str(T1a(j)), 's'))
        hold on
        ylabel('Z Spectra'); ylim([0 1]); ax = gca; ax.XDir = 'reverse';
        h(i) = plot(dwa,Z,'DisplayName',strcat('kb= ',num2str(kb(i)),' Hz'));
        legend(h(1:end),'Location','best')
        subplot(2,2,j+2)
        xline([500 -500]*B0*gamma,'--','Measure range'); xlabel('ppm')
        title(strcat('T1a=',num2str(T1a(j)), 'Hz'))
        hold on; xlim([dwa(end) dwa(1)])
        ylabel('A Spectra');ylim([0 0.8]); ax = gca; ax.XDir = 'reverse';
        plot(dwa(1:domain),A)
    end
end

%% Lithium DEST
% Interaction between electrolyte (Li in organic solvent) and SEI (bound)
% Experiment done on Li-7 spins at 9.4T
gamma = 16.546; % MHz/T; gyromagnetic ratio
B0 = 9.4; % T

% System parameters:
T1a = [10; 100]; T2a = 0.1e-3; % s
T1b = 2.2; T2b = 0.4; % s
M0a = 1; M0b = M0a*5e-3; % arb
kb = [10;50;100;1000;2000]; % s^-1 Select values between 10-1000 Hz
w1 = 750; % Hz (irradiation intensity)
dwa = 2000*B0*gamma:-B0*gamma:-2000*B0*gamma; % Hz (Saturation at given freqs)
db = 0; % offset between wa and wb; dwb = dwa + db

% Calculate spectra and plot results
h = zeros(length(w1),1);
figure
sgtitle('Lithium DEST: Electrolyte + SEI')
for j=1:length(T1a)
    for i=1:length(kb)
        [Z,A,domain] = CEST(T1a(j),T2a,T1b,T2b,kb(i),M0a,M0b,dwa,db,w1);
        subplot(2,2,j)
        xline([500 -500]*B0*gamma,'--'); xlabel('ppm')
        title(strcat('T1a=',num2str(T1a(j)), 's'))
        hold on
        ylabel('Z Spectra'); ylim([0 0.8]); ax = gca; ax.XDir = 'reverse';
        h(i) = plot(dwa,Z,'DisplayName',strcat('kb= ',num2str(kb(i)),' Hz'));
        legend(h(1:end),'Location','best')
        subplot(2,2,j+2)
        xline([500 -500]*B0*gamma,'--','Measure range'); xlabel('ppm')
        title(strcat('T1a=',num2str(T1a(j)), 'Hz'))
        hold on; xlim([dwa(end) dwa(1)])
        ylabel('A Spectra');ylim([0 0.08]); ax = gca; ax.XDir = 'reverse';
        plot(dwa(1:domain),A)
    end
end

%% Sodium CEST
clc;clear;close all;
% Interaction between Na dendrites (free) and SEI (bound)
% Experiment done on Na-23 spins at 9.4T
gamma = 11.262; % MHz/T; gyromagnetic ratio
B0 = 9.4; % T

% System parameters:
T1a = [10; 100]; T2a = 1/(pi*875); % s
T1b = 10e-3; T2b = 4e-3; % s
M0a = 1; M0b = M0a*5e-3; % arb
kb = [10;50;100;500;1000]; % s^-1 Select values between 10-1000 Hz
w1 = 750; % Hz (irradiation intensity)
dwa = 2000*B0*gamma:-B0*gamma:-2000*B0*gamma; % Hz (Saturation at given freqs)
db = -1100*B0*gamma; % offset between wa and wb; dwb = dwa + db

% Calculate spectra and plot results
h = zeros(length(w1),1);
figure
sgtitle('Sodium CEST: Dendrites + SEI')
for j=1:length(T1a)
    for i=1:length(kb)
        [Z,A,domain] = CEST(T1a(j),T2a,T1b,T2b,kb(i),M0a,M0b,dwa,db,w1);
        subplot(2,2,j)
        xline([500 -500]*B0*gamma,'--'); xlabel('ppm')
        title(strcat('T1a=',num2str(T1a(j)), 's'))
        hold on
        ylabel('Z Spectra'); ylim([0 1]); ax = gca; ax.XDir = 'reverse';
        h(i) = plot(dwa,Z,'DisplayName',strcat('kb= ',num2str(kb(i)),' Hz'));
        legend(h(1:end),'Location','best')
        subplot(2,2,j+2)
        xline([500 -500]*B0*gamma,'--','Measure range'); xlabel('ppm')
        title(strcat('T1a=',num2str(T1a(j)), 'Hz'))
        hold on; xlim([dwa(end) dwa(1)])
        ylabel('A Spectra');ylim([0 0.8]); ax = gca; ax.XDir = 'reverse';
        plot(dwa(1:domain),A)
    end
end

%% Soidum DEST
% Interaction between electrolyte (Na in organic solvent) and SEI (bound)
% Experiment done on Na-23 spins at 9.4T
gamma = 11.262; % MHz/T; gyromagnetic ratio
B0 = 9.4; % T

% System parameters:
T1a = [10; 100]; T2a = 0.1e-3; % s
T1b = 3.8e-3; T2b = 3.7e-3; % s
M0a = 1; M0b = M0a*5e-3; % arb
kb = [10;50;100;1000;2000]; % s^-1 Select values between 10-1000 Hz
w1 = 500; % Hz (irradiation intensity)
dwa = 2000*B0*gamma:-B0*gamma:-2000*B0*gamma; % Hz (Saturation at given freqs)
db = 0; % offset between wa and wb; dwb = dwa + db

% Calculate spectra and plot results
h = zeros(length(w1),1);
figure
sgtitle('Sodium DEST: Electrolyte + SEI')
for j=1:length(T1a)
    for i=1:length(kb)
        [Z,A,domain] = CEST(T1a(j),T2a,T1b,T2b,kb(i),M0a,M0b,dwa,db,w1);
        subplot(2,2,j)
        xline([500 -500]*B0*gamma,'--'); xlabel('ppm')
        title(strcat('T1a=',num2str(T1a(j)), 's'))
        hold on
        ylabel('Z Spectra'); ylim([0 0.8]); ax = gca; ax.XDir = 'reverse';
        h(i) = plot(dwa,Z,'DisplayName',strcat('kb= ',num2str(kb(i)),' Hz'));
        legend(h(1:end),'Location','best')
        subplot(2,2,j+2)
        xline([500 -500]*B0*gamma,'--','Measure range'); xlabel('ppm')
        title(strcat('T1a=',num2str(T1a(j)), 'Hz'))
        hold on; xlim([dwa(end) dwa(1)])
        ylabel('A Spectra');ylim([0 0.08]); ax = gca; ax.XDir = 'reverse';
        plot(dwa(1:domain),A)
    end
end
