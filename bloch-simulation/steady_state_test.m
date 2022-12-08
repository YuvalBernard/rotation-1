close all; clc; clear
T1 = [0.1 0.5 1]; % s
T2 = [100e-6 0.5e-3 5e-3]; % s
[T1,T2] = meshgrid(T1,T2);
dw = linspace(0,1000); % Hz 
w1 = linspace(0,1000); % Hz
[DW, W1] = meshgrid(dw,w1);
Mz_ss = zeros(size(DW));
%%%%% Mz_ss %%%%%%
figure
nrows = length(T1); % number of subplot rows
ncols = length(T2); % number of subplot columns
tiledlayout(nrows, ncols)
for i = 1:nrows*ncols
    nexttile()
    for ii = 1:length(w1)
        for jj = 1:length(dw)
            Mz_ss(jj,ii) = blochSS_alt(T1(i),T2(i),W1(jj,ii),DW(jj,ii));
        end
    end
    lvls = 0:0.01:1;
    [c,h] = contourf(W1,DW,Mz_ss,20,'LabelFormat','%0.3f','LevelList',lvls,'LineStyle','None');
    xlabel('w_1 [Hz]')
    ylabel('dw [Hz]')
    title(['Mz^{ss}: T_1 = ',num2str(T1(i)), 's, T_2 = ',num2str(T2(i)),'s'])
end