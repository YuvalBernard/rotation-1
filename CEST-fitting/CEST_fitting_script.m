gamma = 16.546; % MHz/T; gyromagnetic ratio of Li
B0 = 9.4; % T
w0 = gamma*B0; % in MHz

% Read data table.
% cd 'C:\Users\Yuval Bernard\Documents\Weizmann\rotation-1\CEST-fitting\data\derived'
cd 'C:\Users\Yuval Bernard\Documents\Weizmann\rotation-1\CEST-fitting\data\derived'

T = readtable('LP30_dendrotes_CEST_exp_fit.xlsx',...
    'Range','A4:AA55');

xZ = T.Var1 * w0;
Z_500 = T.Var20;
Z_1000 = T.Var22;
Z_1500 = T.Var24;
Z_2000 = T.Var26;

Z = [Z_500,Z_1000,Z_1500,Z_2000];

% Fit from cmdstanr::pml
dendrite = struct('T1',1/8,'T2',1/393);
SEI = struct('T1',1/10.9,'T2',1/(30305),'k',299,'f',0.0121,'dw',-260*w0);
sys = struct('offsets',xZ,'tp',0.2);
w1 = [500,1000,1500,2000];

% Simulate data sets from fit.
figure; ylabel('Z Spectra'); xlabel('offset (ppm)'); ylim([0 1]); set(gca,'XDir','reverse');
h = zeros(size(w1)); % Placeholder array
C = lines(length(w1)); % Cell array of colros.
for i=1:length(w1)
        sys.w1 = w1(i);
        hold on
        Z_sim = CEST_multipool(sys,dendrite,SEI);
        plot(xZ/w0,Z(:,i),'.','Color',C(i,:))
        h(i) = plot(xZ/w0,Z_sim,'Color',C(i,:),'DisplayName',['w1=',num2str(sys.w1),' Hz']);
end
legend(h(1:end),'Location','Southeast');
