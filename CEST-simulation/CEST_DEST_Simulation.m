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

figure;
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

f = figure('Visible','off');
for i = 1:length(w1)
    sys.w1 = w1(i);    
    Z = CEST_multipool(sys,dendrite,SEI);
    hold on
    ylabel('Z Spectra'); ylim([0 1]); ax = gca; ax.XDir = 'reverse';
    plot(sys.offsets/w0,Z);
end
set(f,'Visible','on')
% f.Position = [721,311,304.5,331.5];
%% Lithium CEST
clc;clearvars t sys SEI dendrite ;
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
for j=1:length(SEI_T1)
    SEI.T1 = SEI_T1(j);
    nexttile
    for i=1:length(SEI_k)
        SEI.k = SEI_k(i);
        title(strcat('T_1 [SEI]=',num2str(SEI.T1), 's'))
        hold on
        ylabel('Z Spectra'); ylim([0 1]); set(gca,'XDir','reverse')
        [Z,MTR_asym] = CEST_multipool(sys,dendrite,SEI);
        h(i) = plot(sys.offsets/w0,Z,'DisplayName',['k=',num2str(SEI.k),' Hz']);
        plot(sys.offsets/w0,MTR_asym,'--')
    end
    legend(h(1:end),'Location','Southeast');
    
end
%% Lithium DEST
% Interaction between electrolyte (Li in organic solvent) and SEI (bound)
% Experiment done on Li-7 spins at 9.4T
gamma = 16.546; % MHz/T; gyromagnetic ratio
B0 = 9.4; % T
w0 = gamma*B0;

% System parameters:
electrolyte = struct('T1',2.2,'T2',0.4);
SEI = struct('T1',[10,100],'T2',0.1e-3,'f',0.1,'dw',0);
sys = struct('w1',750,'tp',5,'offsets',linspace(-500,500,5000)*w0);
SEI_T1 = [10;100]; SEI_k = [10;50;100;500;1000];
% Calculate spectra and plot results
f = figure('Visible','off');
t = tiledlayout(1,length(SEI_T1));
t.TileSpacing = 'compact'; t.Padding = 'compact';
title(t,'Lithium DEST: Electrolyte + SEI',...
    ['measure at +-500ppm, \omega_1 = ',num2str(sys.w1),' Hz'])
xlabel(t,'\Delta\omega [ppm]')
h = zeros(size(SEI_k));
for j=1:length(SEI_T1)
    SEI.T1 = SEI_T1(j);
    nexttile
    for i=1:length(SEI_k)
        SEI.k = SEI_k(i);
        title(strcat('T_1 [SEI]=',num2str(SEI.T1), 's'))
        hold on
        ylabel('Z Spectra'); ylim([0 1]); set(gca,'XDir','reverse')
        [Z,MTR] = CEST_multipool(sys,electrolyte,SEI);
        h(i) = plot(sys.offsets/w0,Z,'DisplayName',['k=',num2str(SEI.k),' Hz']);
        plot(sys.offsets/w0,MTR,'--')
    end
    legend(h(1:end),'Location','Southeast');
end
set(f,'Visible','on')
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
f = figure('Visible','off');
t = tiledlayout(1,length(SEI_T1));
t.TileSpacing = 'compact'; t.Padding = 'compact';
title(t,'Sodium CEST: Dendrites + SEI',...
    ['measure at +-1500ppm, \omega_1 = ',num2str(sys.w1),' Hz'])
xlabel(t,'\Delta\omega [ppm]')
h = zeros(size(SEI_k));
for j=1:length(SEI_T1)
    SEI.T1 = SEI_T1(j);
    nexttile
    for i=1:length(SEI_k)
        SEI.k = SEI_k(i);
        title(strcat('T_1 [SEI]=',num2str(SEI.T1), 's'))
        hold on
        ylabel('Z Spectra'); ylim([0 1]); set(gca,'XDir','reverse')
        h(i) = plot(sys.offsets/w0,CEST_multipool(sys,dendrite,SEI),'DisplayName',['k=',num2str(SEI.k),' Hz']);
    end
    legend(h(1:end),'Location','Southeast');
end
set(f,'Visible','on')
%% Soidum DEST
% Interaction between electrolyte (Na in organic solvent) and SEI (bound)
% Experiment done on Na-23 spins at 9.4T
gamma = 11.262; % MHz/T; gyromagnetic ratio
B0 = 9.4; % T
w0 = gamma*B0;

% System parameters:
electrolyte = struct('T1',3.8e-3,'T2',3.7e-3);
SEI = struct('T1',[10,100],'T2',0.1e-3,'f',0.02,'dw',0);
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
for j=1:length(SEI_T1)
    SEI.T1 = SEI_T1(j);
    nexttile
    for i=1:length(SEI_k)
        SEI.k = SEI_k(i);
        title(strcat('T_1 [SEI]=',num2str(SEI.T1), 's'))
        hold on
        ylabel('Z Spectra'); ylim([0 1]); set(gca,'XDir','reverse')
        h(i) = plot(sys.offsets/w0,CEST_multipool(sys,electrolyte,SEI),'DisplayName',['k=',num2str(SEI.k),' Hz']);
    end
    legend(h(1:end),'Location','Southeast');
end


%% Limit Testing: Li CEST
% Li CEST: What are the limits of the exchanges rates that can be measured
% Li dendrites/SEI (fraction of 1:0.02)

% Given domains: SEI_k in [0.01 1000] Hz and sys_w1 in [250 2000] Hz.

% Method: System is called analyzable if direct solving of Bloch-McConnell
% equations yields observable, significant peaks in asymmetry spectra.
% Significant is rather subjective. We limit peaks to be at least 0.05 tall

% Create mesh grid of different kb and w1 values. Check analyzability
% criterion for all combinations in grid and create map of peak height.
clc;clear;close all;

% System parameters
gamma = 16.546; % MHz/T; gyromagnetic ratio
B0 = 9.4; % T
w0 = gamma*B0; % MHz

dendrite = struct('T1',150e-3,'T2',0.5e-3);
SEI = struct('T1',10,'T2',1/(35e3),'f',0.02,'dw',-260*w0);
sys = struct('tp',1,'offsets',linspace(-750,750,1000)*w0);

% Create grid: 
[SEI_k,sys_w1] = meshgrid(linspace(0.01,1000,25),linspace(250,2000,25));

% Initialize peaks array
pks = zeros(size(SEI_k));
f = figure('Visible','off');
for jj = 1:size(sys_w1,2)
    for ii = 1:size(SEI_k,1)
        % Calculate spectra
        SEI.k = SEI_k(ii,jj); sys.w1 = sys_w1(ii,jj);
        [~,MTR] = CEST_multipool(sys,dendrite,SEI);
        % Find peaks in MTR spectrum
        % Robust method for more than one peak
%         pks(ii,jj) = ~isempty(findpeaks(A),'MinPeakHeight',0.05));
        % Alternative method, less robust per se
        pks(ii,jj) = max(MTR)*(max(MTR) >= 0.05);
    end
end

% Plot results
contourf(SEI_k,sys_w1,pks,'LineStyle','None')
xlabel('k_{SEI} [Hz]'); ylabel('\omega_1 [Hz]'); title("Li CEST: Peaks height in MTR spectra"); colorbar
set(f,'Visible','on')

%% Limit Testing: Na CEST
%  CEST: What are the limits of the exchanges rates that can be measured
% Na dendrites/SEI (fraction of 1:0.02)

% Given domains: kb in [0.01 1000] Hz and w1 in [250 2000] Hz.

% Create mesh grid of different kb and w1 values. Check analyzability
% criterion for all combinations in grid and create map of peak height.
clc;clear;close all;

% System parameters
gamma = 16.546; % MHz/T; gyromagnetic ratio
B0 = 9.4; % T
w0 = gamma*B0; % MHz

dendrite = struct('T1',10e-3,'T2',4e-3);
SEI = struct('T1',10,'T2',1/(875e3),'f',0.02,'dw',-1100*w0);
sys = struct('tp',1,'offsets',linspace(-1500,1500,1000)*w0);

% Create grid: 
[SEI_k,sys_w1] = meshgrid(linspace(0.01,1000,25),linspace(250,2000,25));

% Initialize peaks array
pks = zeros(size(SEI_k));

for jj = 1:size(sys_w1,2)
    for ii = 1:size(SEI_k,1)
        % Calculate spectra
        SEI.k = SEI_k(ii,jj); sys.w1 = sys_w1(ii,jj);
        [~,MTR] = CEST_multipool(sys,dendrite,SEI);
        % Find peaks in MTR spectrum
        % Robust method for more than one peak
%         pks(ii,jj) = ~isempty(findpeaks(A),'MinPeakHeight',0.05));
        % Alternative method, less robust per se
        pks(ii,jj) = max(MTR)*(max(MTR) >= 0.05);
    end
end

% Plot results
contourf(SEI_k,sys_w1,pks,'LineStyle','None')
xlabel('k_{SEI} [Hz]'); ylabel('\omega_1 [Hz]'); title("Na CEST: Peaks height in MTR spectra"); colorbar

%% Limit Testing: Li DEST
% Li DEST: What are the limits of the exchanges rates that can be measured
% Li electrolyte/SEI.

% Given domains: SEI_k in [0.01 1000] Hz and SEI_f in [0.01 0.1].
% Calculate for sys_w1 = [250, 500, 2000] Hz.

% Method: System is called analyzable if direct solving of Bloch-McConnell
% equations yields observable, significant peaks in asymmetry spectra.
% Significant is rather subjective. We limit peaks to be at least 0.05 tall

% Create mesh grid of different kb and w1 values. Check analyzability
% criterion for all combinations in grid and create map of peak height.
clc;clear;close all;

% System parameters
gamma = 16.546; % MHz/T; gyromagnetic ratio
B0 = 9.4; % T
w0 = gamma*B0; % MHz

electrolyte = struct('T1',2.2,'T2',0.4);
SEI = struct('T1',10,'T2',0.1e-3,'dw',0);
sys = struct('tp',10,'offsets',linspace(-750,750,1000)*w0);

% Create grid: 
[SEI_k,SEI_f] = meshgrid(linspace(0.01,1000,10),linspace(0.01,0.1,10));
% Create sys_w1 array
sys_w1 = [250;500;2000];

% Initialize peaks array
pks = zeros(size(SEI_k));
f = figure('Visible','off');

t = tiledlayout(1,length(sys_w1));
t.TileSpacing = 'compact'; t.Padding = 'compact';
title(t,'Lithium DEST: Electrolyte + SEI',...
    'measure at +-750ppm, Varying \omega_1')
xlabel(t,'\Delta\omega [ppm]')
for kk = 1:length(sys_w1)
    sys.w1 = sys_w1(kk);
    for jj = 1:size(SEI_f,2)
        for ii = 1:size(SEI_k,1)
            % Calculate spectra
            SEI.k = SEI_k(ii,jj); SEI.f = SEI_f(ii,jj);
            [~,MTR] = CEST_multipool(sys,electrolyte,SEI);
            % Find peaks in MTR spectrum
            pks(ii,jj) = max(MTR)*(max(MTR) >= 0.05);
        end
    end
    nexttile
    % Plot results
    contourf(SEI_k,SEI_f,pks,'LineStyle','None')
    xlabel('k_{SEI} [Hz]'); ylabel('f_{SEI}'); title(['\omega_1 = ',num2str(sys.w1),' Hz']); colorbar
end
set(f,'Visible','on')