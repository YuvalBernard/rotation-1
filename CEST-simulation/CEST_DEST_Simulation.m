% CEST and DEST simulations (Slides from 19.12.22)


%% Lithium CEST: Reconstruct experimental results:  DOI: 10.1021/jacs.2c02494 
% LP30/FEC @ 298K
clc;clearvars gamma B0 w0 dendrite SEI sys w1; close all

gamma = 16.546; % MHz/T; gyromagnetic ratio
B0 = 9.4; % T
w0 = gamma*B0; % in MHz

dendrite = struct('T1',1/7,'T2',1/1456);
SEI = struct('T1',100,'T2',1/(28e3),'k',91,'f',0.02,'dw',-260*w0);
sys = struct('offsets',(500:-1:-500)*w0,'tp',0.2);
w1 = [500; 1000; 1500; 2000];

f = figure;
for i = 1:length(w1)
    sys.w1 = w1(i);    
    Z = CEST_multipool(sys,dendrite,SEI);
    hold on
    ylabel('Z Spectra'); ylim([0 1]); ax = gca; ax.XDir = 'reverse';
    plot(sys.offsets/w0,Z);
end

%%
% LP30 @ 323K
clc;clearvars gamma B0 w0 dendrite SEI sys w1; close all

gamma = 16.546; % MHz/T; gyromagnetic ratio
B0 = 9.4; % T
w0 = gamma*B0; % in MHz

dendrite = struct('T1',1/8,'T2',1/393);
SEI = struct('T1',100,'T2',1/(28e3),'k',285,'f',0.02,'dw',-260*w0);
sys = struct('offsets',(500:-1:-500)*w0,'tp',0.2);
w1 = [500; 1000; 1500; 2000];

f = figure;
for i = 1:length(w1)
    sys.w1 = w1(i);    
    Z = CEST_multipool(sys,dendrite,SEI);
    hold on
    ylabel('Z Spectra'); ylim([0 1]); ax = gca; ax.XDir = 'reverse';
    plot(sys.offsets/w0,Z);
end
% f.Position = [721,311,304.5,331.5];
%% Lithium CEST
clc;clear;close all;
% Interaction between Li dendrites (free) and SEI (bound)
% Experiment done on Li-7 spins at 9.4T
gamma = 16.546; % MHz/T; gyromagnetic ratio
B0 = 9.4; % T
w0 = gamma*B0;

% System parameters:
dendrite = struct('T1',150e-3,'T2',0.5e-3);
SEI = struct('T1',[10,100],'T2',1/(35e3),'f',0.02,'dw',-260*w0);
sys = struct('w1',750,'tp',0.2,'offsets',linspace(-500,500,1000)*w0);
SEI_T1 = [10;100]; SEI_k = [10;50;100;500;1000];
% Calculate spectra and plot results
figure
t = tiledlayout(1,length(SEI_T1));
t.TileSpacing = 'compact'; t.Padding = 'compact';
title(t,'Lithium CEST: Dendrites + SEI',...
    ['measure at +-500ppm, \omega_1 = ',num2str(sys.w1),' Hz'])
xlabel(t,'\Delta\omega [ppm]')
h = zeros(size(SEI_k));
nexttile
for j=1:length(SEI_T1)
    SEI.T1 = SEI_T1(j);
    for i=1:length(SEI_k)
        SEI.k = SEI_k(i);
        title(strcat('T_1 [SEI]=',num2str(SEI.T1), 's'))
        hold on
        ylabel('Z Spectra'); ylim([0 1]); set(gca,'XDir','reverse')
        h(i) = plot(sys.offsets/w0,CEST_multipool(sys,dendrite,SEI),'DisplayName',['k=',num2str(SEI.k),' Hz']);
    end
    legend(h(1:end),'Location','Southeast');
    nexttile
end
%% Lithium DEST
% Interaction between electrolyte (Li in organic solvent) and SEI (bound)
% Experiment done on Li-7 spins at 9.4T
gamma = 16.546; % MHz/T; gyromagnetic ratio
B0 = 9.4; % T
w0 = gamma*B0;

% System parameters:
electrolyte = struct('T1',2.2,'T2',0.4);
SEI = struct('T1',[10,100],'T2',0.1e-3,'f',0.02,'dw',0);
sys = struct('w1',750,'tp',5,'offsets',linspace(-500,500,5000)*w0);
SEI_T1 = [10;100]; SEI_k = [10;50;100;500;1000];
% Calculate spectra and plot results
figure
t = tiledlayout(1,length(SEI_T1));
t.TileSpacing = 'compact'; t.Padding = 'compact';
title(t,'Lithium DEST: Electrolyte + SEI',...
    ['measure at +-500ppm, \omega_1 = ',num2str(sys.w1),' Hz'])
xlabel(t,'\Delta\omega [ppm]')
h = zeros(size(SEI_k));
nexttile
for j=1:length(SEI_T1)
    SEI.T1 = SEI_T1(j);
    for i=1:length(SEI_k)
        SEI.k = SEI_k(i);
        title(strcat('T_1 [SEI]=',num2str(SEI.T1), 's'))
        hold on
        ylabel('Z Spectra'); ylim([0 1]); set(gca,'XDir','reverse')
        h(i) = plot(sys.offsets/w0,CEST_multipool(sys,electrolyte,SEI),'DisplayName',['k=',num2str(SEI.k),' Hz']);
    end
    legend(h(1:end),'Location','Southeast');
    nexttile
end
%% Sodium CEST
clc;clear;close all;
% Interaction between Na dendrites (free) and SEI (bound)
% Experiment done on Na-23 spins at 9.4T
gamma = 11.262; % MHz/T; gyromagnetic ratio
B0 = 9.4; % T
w0 = gamma*B0;

% System parameters:
dendrite = struct('T1',10e-3,'T2',4e-3);
SEI = struct('T1',[10,100],'T2',1/(875e3),'f',0.02,'dw',-1100*w0);
sys = struct('w1',750,'tp',5,'offsets',linspace(-1500,1500,5000)*w0);
SEI_T1 = [10;100]; SEI_k = [10;50;100;500;1000];
% Calculate spectra and plot results
figure
t = tiledlayout(1,length(SEI_T1));
t.TileSpacing = 'compact'; t.Padding = 'compact';
title(t,'Sodium CEST: Dendrites + SEI',...
    ['measure at +-1500ppm, \omega_1 = ',num2str(sys.w1),' Hz'])
xlabel(t,'\Delta\omega [ppm]')
h = zeros(size(SEI_k));
nexttile
for j=1:length(SEI_T1)
    SEI.T1 = SEI_T1(j);
    for i=1:length(SEI_k)
        SEI.k = SEI_k(i);
        title(strcat('T_1 [SEI]=',num2str(SEI.T1), 's'))
        hold on
        ylabel('Z Spectra'); ylim([0 1]); set(gca,'XDir','reverse')
        h(i) = plot(sys.offsets/w0,CEST_multipool(sys,dendrite,SEI),'DisplayName',['k=',num2str(SEI.k),' Hz']);
    end
    legend(h(1:end),'Location','Southeast');
    nexttile
end

%% Soidum DEST
% Interaction between electrolyte (Na in organic solvent) and SEI (bound)
% Experiment done on Na-23 spins at 9.4T
gamma = 11.262; % MHz/T; gyromagnetic ratio
B0 = 9.4; % T
w0 = gamma*B0;

% System parameters:
electrolyte = struct('T1',3.8e-3,'T2',3.7e-3);
SEI = struct('T1',[10,100],'T2',0.1,'f',0.02,'dw',0);
sys = struct('w1',750,'tp',5,'offsets',linspace(-1500,1500,5000)*w0);
SEI_T1 = [10;100]; SEI_k = [10;50;100;500;1000];
% Calculate spectra and plot results
figure
t = tiledlayout(1,length(SEI_T1));
t.TileSpacing = 'compact'; t.Padding = 'compact';
title(t,'Sodium DEST: Electrolyte + SEI',...
    ['measure at +-1500ppm, \omega_1 = ',num2str(sys.w1),' Hz'])
xlabel(t,'\Delta\omega [ppm]')
h = zeros(size(SEI_k));
nexttile
for j=1:length(SEI_T1)
    SEI.T1 = SEI_T1(j);
    for i=1:length(SEI_k)
        SEI.k = SEI_k(i);
        title(strcat('T_1 [SEI]=',num2str(SEI.T1), 's'))
        hold on
        ylabel('Z Spectra'); ylim([0 1]); set(gca,'XDir','reverse')
        h(i) = plot(sys.offsets/w0,CEST_multipool(sys,electrolyte,SEI),'DisplayName',['k=',num2str(SEI.k),' Hz']);
    end
    legend(h(1:end),'Location','Southeast');
    nexttile
end


%% Limit Testing: Li CEST
% Li CEST: What are the limits of the exchanges rates that can be measured
% Li dendrites/SEI (fraction of 1:0.02)

% Given domains: kb in [0.01 1000] Hz and w1 in [250 2000] Hz.

% Method: System is called analyzable if direct solving of Bloch-McConnell
% equations yields observable, significant peaks in asymmetry spectra.
% Significant is rather subjective. We limit peaks to be at least 0.05 tall

% Create mesh grid of different kb and w1 values. Check analyzability
% criteria for all combinations in grid and create map of peak height.
clc;clear;close all;

% System parameters
gamma = 16.546; % MHz/T; gyromagnetic ratio
B0 = 9.4; % T
w0 = gamma*B0; % MHz

T1a = 10; T2a = 1/(pi*35); % s
T1b = 150e-3; T2b = 0.5e-3; % s
M0a = 1; M0b = 0.02; % arb
dwa = 1000*w0:-w0:-1000*w0; % Hz (Saturation at given freqs)
db = -260*w0; % offset between wa and wb; dwb = dwa + db

% Create grid: 
[kb,w1] = meshgrid(linspace(0.01,1000,25),linspace(250,2000,25));

% Initialize peaks array
pks = zeros(size(kb));

for jj = 1:length(w1)
    for ii = 1:length(kb)
        % Calculate spectra
        [~,A,domain] = CEST(T1a,T2a,T1b,T2b,kb(ii,jj),M0a,M0b,dwa,db,w1(ii,jj),1);
        % Find peaks in Asymmetry spectrum
        % Robust method for more than one peak
%         pks(ii,jj) = ~isempty(findpeaks(A),'MinPeakHeight',0.05));
        % Alternative method, less robust per se
        pks(ii,jj) = max(A)*(max(A) >= 0.05);
    end
end

% Plot results
lvls = 0:0.01:1;
contourf(kb,w1,pks,'LabelFormat','%0.3f','LevelList',lvls,'LineStyle','None')
xlabel('k_b [Hz]'); ylabel('\omega_1 [Hz]'); title("Li CEST: Peaks' height in Asymmetry spectra"); colorbar
% lvls = 0:0.01:1;
% contourf(kb,w1,pks,20,'LabelFormat','%0.3f','LevelList',lvls,'LineStyle','None');

% % Inspect results for single w1,kb
% w1 = 750; kb = 10;
% [Z,A,domain] = CEST(T1a,T2a,T1b,T2b,kb,M0a,M0b,dwa,db,w1);
% t = tiledlayout(2,1);
% title(t,['T_{1a} = ',num2str(T1a),'s,' ...
%     ' \omega_1 = ',num2str(w1),'Hz, k_b = ',num2str(kb),'Hz'])
% ax1 = nexttile;
% plot(dwa,Z); ax = gca; ax.XDir = 'reverse'; title('Z spectrum')
% ax2 = nexttile; hold on; linkaxes([ax1,ax2],'x');
% plot(dwa(1:domain),A); ax = gca; ax.XDir = 'reverse';  title('A spectrum')
% 
% % Find peaks in Asymmetry spectrum
% pks = findpeaks(A,fliplr(dwa(1:domain)),'MinPeakHeight',0.05);
% plot(dwa(pks == A),pks,'o'); ax = gca; ax.XDir = 'reverse';
