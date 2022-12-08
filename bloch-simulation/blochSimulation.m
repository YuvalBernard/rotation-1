%% Plot Mx(t),My(t),Mz(t); if you want to look at a specific case.
clear; clc; close all
% Define constants: default values
T1 = 1; % s
T2 = 1e-3; %s
w1 = 150; % Hz; gamma*H1
dw = 2500; % Hz; gamma*H0-w
t0 = 0; tmax = 7; q = 10000; t = t0:(tmax-t0)/(q-1):tmax;
[Mz,Mx,My] = bloch(T1,T2,w1,dw,t0,tmax,q);
% [t_ss,Mz_ss] = blochSS_alt(T1,T2,w1,dw,t0,tmax,q,Mz);
[t_ss,Mz_ss] = blochSS(t,Mz);
figure; plot(t,Mz,'b-','LineWidth',2); hold on; plot(t_ss,Mz_ss,'r.','MarkerSize',24)
figure; plot(t,Mx,t,My,t,Mz); legend('Mx','My','Mz')

%% Compare numerical and analytical solutions
R1 = 1/T1; % s^-1
R2 = 1/T2; % s^-1
Mz0 = 1; % initial z-magnetization
Anal_Mz_ss = Mz0*(R1*(R2^2+dw^2))/(R1*(R2^2+dw^2)+w1^2*R2);

%% Change T1
% close all; clc; clear
T2 = 2e-3;
w1 = 150;
dw = 2500;
t0 = 0; tmax = 7; q = 10000; t = t0:(tmax-t0)/(q-1):tmax;
T1 = [0.1,0.5,1,2,10];
colors = lines(length(T1));
h = zeros(size(T1));
figure(); hold on; grid on; xlabel('t (s)'); ylabel('Mz'); 
title(['T_2 = ',num2str(T2),'s, w_1 = ', num2str(w1),'Hz, dw = ',num2str(dw),'Hz'])
for i = 1:length(T1)
    Mz = bloch(T1(i),T2,w1,dw,t0,tmax,q);
%     [t_ss,Mz_ss] = blochSS(t,Mz);
    [t_ss,Mz_ss] = blochSS_alt(T1(i),T2,w1,dw,t0,tmax,q,Mz);
    h(i) = plot(t,Mz,'color',colors(i,:),'LineWidth',2,'DisplayName',strcat('T1= ',num2str(T1(i))));
    plot(t_ss,Mz_ss,'.','color',colors(i,:),'MarkerSize',24)
end
legend(h(1:end))

%%%%%% Conclusion: Increaing T1 decreases Mz_ss and increases t_ss. %%%%%%
%% Change T2
close all; clc; clear
T1 = 1;
w1 = 150;
dw = 2500;
t0 = 0; tmax = 7; q = 10000; t = t0:(tmax-t0)/(q-1):tmax;
T2 = [20e-3, 5e-3, 1e-3, 500e-6, 100e-6, 10e-6];
colors = lines(length(T2));
h = zeros(size(T2));
figure(); hold on; grid on; xlabel('t (s)'); ylabel('Mz'); 
title(['T_1 = ',num2str(T1),'s, w_1 = ', num2str(w1),'Hz, dw = ',num2str(dw),'Hz'])
for i = 1:length(T2)
    Mz = bloch(T1,T2(i),w1,dw,t0,tmax,q);
    [t_ss,Mz_ss] = blochSS(t,Mz);
    h(i) = plot(t,Mz,'color',colors(i,:),'LineWidth',2,'DisplayName',strcat('T2= ',num2str(T2(i))));
    plot(t_ss,Mz_ss,'.','color',colors(i,:),'MarkerSize',24)
end
legend(h(1:end))

%%%%%% Conclusion: Increasing T2 increases Mz_ss and t_ss. %%%%%%
%% Change w1
close all; clc; clear
T1 = 1;
T2 = 5e-3;
dw = 2500;
t0 = 0; tmax = 7; q = 10000; t = t0:(tmax-t0)/(q-1):tmax;
w1 = [50,150,300,750,1500];
colors = lines(length(w1));
h = zeros(size(w1));
figure(); hold on; grid on; xlabel('t (s)'); ylabel('Mz'); 
title(['T_2 = ',num2str(T2),'s, T_1 = ', num2str(T1),'s, dw = ',num2str(dw),'Hz'])
for i = 1:length(w1)
    Mz = bloch(T1,T2,w1(i),dw,t0,tmax,q);
    [t_ss,Mz_ss] = blochSS(t,Mz);
    h(i) = plot(t,Mz,'color',colors(i,:),'LineWidth',2,'DisplayName',strcat('w_1= ',num2str(w1(i))));
    plot(t_ss,Mz_ss,'.','color',colors(i,:),'MarkerSize',24)
end
legend(h(1:end))

%%%%%% Conclusion: Increasing w1 decreases Mz_ss and t_ss. %%%%%%
%% Change dw
close all; clc; clear
T1 = 1;
T2 = 5e-3;
w1 = 150;
t0 = 0; tmax = 7; q = 10000; t = t0:(tmax-t0)/(q-1):tmax;
dw = [0,250,1000,2500,5000];
colors = lines(length(dw));
h = zeros(size(dw));
figure(); hold on; grid on; xlabel('t (s)'); ylabel('Mz'); 
title(['T_2 = ',num2str(T2),'s, T_1 = ', num2str(T1),'s, w_1 = ',num2str(w1),'Hz'])
for i = 1:length(dw)
    Mz = bloch(T1,T2,w1,dw(i),t0,tmax,q);
    [t_ss,Mz_ss] = blochSS(t,Mz);
    h(i) = plot(t,Mz,'color',colors(i,:),'LineWidth',2,'DisplayName',strcat('dw= ',num2str(dw(i))));
    plot(t_ss,Mz_ss,'.','color',colors(i,:),'MarkerSize',24)

end
legend(h(1:end))

%%%%%% Conclusion: Increasing dw increases Mz_ss and t_ss. %%%%%%
%% Mz_ss and t_ss dependence on w1,dw (for specific T1,T2)
tic
close all; clc; clear
T1 = [0.1 0.5 1]; % s
T2 = [100e-6 0.5e-3 5e-3]; % s
[T1,T2] = meshgrid(T1,T2);
t0 = 0; tmax = 7; q = 100; t = t0:(tmax-t0)/(q-1):tmax; % s
dw = linspace(0,1000,100); % Hz 
w1 = linspace(0,1000,100); % Hz
[DW, W1] = meshgrid(dw,w1);
Mz_ss = zeros(size(DW));
%%%%% Mz_ss %%%%%%
figure
nrows = 3; % number of subplot rows
ncols = 3; % number of subplot columns
tiledlayout(nrows, ncols)
for i = 1:nrows*ncols
    nexttile()
    for ii = 1:length(w1)
        for jj = 1:length(dw)
            [~,Mz_ss(jj,ii)] = blochSS_alt(T1(i),T2(i),W1(jj,ii),DW(jj,ii),t0,tmax,q,0);
%             Mz_ss(jj,ii) = blochSS(t,bloch(T1(i),T2(i),W1(jj,ii),DW(jj,ii),t0,tmax,q));
        end
    end
    lvls = 0:0.01:1;
    [c,h] = contourf(W1,DW,Mz_ss,20,'LabelFormat','%0.3f','LevelList',lvls,'LineStyle','None');
    xlabel('w_1 [Hz]')
    ylabel('dw [Hz]')
    title(['Mz^{ss}: T_1 = ',num2str(T1(i)), 's, T_2 = ',num2str(T2(i)),'s'])
end
%%%% t_ss %%%%%%
figure
nrows = 3; % number of subplot rows
ncols = 3; % number of subplot columns
tiledlayout(nrows, ncols)
t_ss = zeros(size(DW));
for i = 1:nrows*ncols
    nexttile()
    for ii = 1:length(w1)
        for jj = 1:length(dw)
            Mz = bloch(T1(i),T2(i),W1(jj,ii),DW(jj,ii),t0,tmax,q);
%             t_ss(jj,ii) = blochSS_alt(T1(i),T2(i),W1(jj,ii),DW(jj,ii),t0,tmax,q,Mz);
            t_ss(jj,ii) = blochSS(t,Mz);
        end
    end
    [c,h] = contourf(W1,DW,t_ss,20,'LabelFormat','%0.3f','LevelList',lvls,'LineStyle','None');
    xlabel('w_1 [Hz]')
    ylabel('dw [Hz]')
    title(['t^{ss}: T_1 = ',num2str(T1(i)), 's, T_2 = ',num2str(T2(i)),'s'])
end
toc