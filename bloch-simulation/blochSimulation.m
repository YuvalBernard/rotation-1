%% Plot Mx(t),My(t),Mz(t); if you want to look at a specific case.
clear; clc; close all
% Define constants: default values
T1 = 1; % s
T2 = 5e-3; %s
w1 = 750; % Hz; gamma*H1
dw = 2500; % Hz; gamma*H0-w
t = linspace(0,7,10000)';
[Mz,Mx,My] = bloch(t,T1,T2,w1,dw);
[t_ss,Mz_ss] = blochSS(t,Mz);
figure; plot(t,Mz,'b-','LineWidth',2); hold on; plot(t_ss,Mz_ss,'r.','MarkerSize',24)
figure; plot(t,Mx,t,My,t,Mz); legend('Mx','My','Mz')

%% Compare numerical and analytical solutions
R1 = 1/T1; % s^-1
R2 = 1/T2; % s^-1
Mz0 = 1; % initial z-magnetization
Anal_Mz_ss = Mz0*(R1*(R2^2+dw^2))/(R1*(R2^2+dw^2)+w1^2*R2);

%% Change T1
close all; clc; clear
T2 = 2e-3;
w1 = 150;
dw = 2500;
t = linspace(0,7,10000)';
T1 = [0.1,0.5,1,2,10];
colors = lines(length(T1));
h = zeros(size(T1));
figure(); hold on; grid on; xlabel('t (s)'); ylabel('Mz'); 
title(['T_2 = ',num2str(T2),'s, w_1 = ', num2str(w1),'Hz, dw = ',num2str(dw),'Hz'])
for i = 1:length(T1)
    Mz = bloch(t,T1(i),T2,w1,dw);
    [t_ss,Mz_ss] = blochSS(t,Mz);
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
t = linspace(0,7,10000)';
T2 = [20e-3, 5e-3, 1e-3, 500e-6, 100e-6, 10e-6];
colors = lines(length(T2));
h = zeros(size(T2));
figure(); hold on; grid on; xlabel('t (s)'); ylabel('Mz'); 
title(['T_1 = ',num2str(T1),'s, w_1 = ', num2str(w1),'Hz, dw = ',num2str(dw),'Hz'])
for i = 1:length(T2)
    Mz = bloch(t,T1,T2(i),w1,dw);
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
t = linspace(0,7,10000)';
w1 = [50,150,300,750,1500];
colors = lines(length(w1));
h = zeros(size(w1));
figure(); hold on; grid on; xlabel('t (s)'); ylabel('Mz'); 
title(['T_2 = ',num2str(T2),'s, T_1 = ', num2str(T1),'s, dw = ',num2str(dw),'Hz'])
for i = 1:length(w1)
    Mz = bloch(t,T1,T2,w1(i),dw);
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
t = linspace(0,7,10000)';
dw = [0,250,1000,2500,5000];
colors = lines(length(dw));
h = zeros(size(dw));
figure(); hold on; grid on; xlabel('t (s)'); ylabel('Mz'); 
title(['T_2 = ',num2str(T2),'s, T_1 = ', num2str(T1),'s, w_1 = ',num2str(w1),'Hz'])
for i = 1:length(dw)
    Mz = bloch(t,T1,T2,w1,dw(i));
    [t_ss,Mz_ss] = blochSS(t,Mz);
    h(i) = plot(t,Mz,'color',colors(i,:),'LineWidth',2,'DisplayName',strcat('dw= ',num2str(dw(i))));
    plot(t_ss,Mz_ss,'.','color',colors(i,:),'MarkerSize',24)

end
legend(h(1:end))

%%%%%% Conclusion: Increasing dw increases Mz_ss and t_ss. %%%%%%
%% Mz_ss and t_ss dependence on w1,dw (for specific T1,T2)
close all; clc; clear
T1 = 1;
T2 = 20e-3;
t = linspace(0,7,100)';
dw = linspace(0,2500); 
w1 = linspace(0,150); 
[DW, W1] = meshgrid(dw,w1);
[t_ss,Mz_ss] = arrayfun(@(w1,dw) blochSS(t,bloch(t,T1,T2,w1,dw)), W1, DW);
%%%%% Mz_ss %%%%%%

% figure
% meshc(W1,DW,Mz_ss)
% xlabel('w_1 [Hz]')
% ylabel('dw [Hz]')
% zlabel('Mz^{ss}')
% title(['Mz^{ss}: ''T_1 = ',num2str(T1), 's, T_2 = ',num2str(T2),'s'])

figure
contourf(W1,DW,Mz_ss,'ShowText',true,'LabelFormat','%0.3f','LevelList', 0:max(Mz_ss,[],'all'))
xlabel('w_1 [Hz]')
ylabel('dw [Hz]')
zlabel('Mz^{ss}')
title(['Mz^{ss}: ''T_1 = ',num2str(T1), 's, T_2 = ',num2str(T2),'s'])

%%%%% t_ss %%%%%%

% figure
% meshc(W1,DW,t_ss)
% xlabel('w_1 [Hz]')
% ylabel('dw [Hz]')
% zlabel('t^{ss} [s]')
% title(['t^{ss}: ''T_1 = ',num2str(T1), 's, T_2 = ',num2str(T2),'s'])

figure
contourf(W1,DW,t_ss,"ShowText",true,"LabelFormat","%0.3f",'LevelList', 0:max(t_ss,[],'all'))
xlabel('w_1 [Hz]')
ylabel('dw [Hz]')
zlabel('t^{ss} [s]')
title(['t^{ss}: ''T_1 = ',num2str(T1), 's, T_2 = ',num2str(T2),'s'])

% [opt_w1, opt_dw, opt_t_ss, opt_Mz_ss, Summary] = blochOptimalSaturationPulse(t,T1,T2,w1,dw);
% disp(Summary);
