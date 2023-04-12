
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% this file is a skeleton of a script that you may use as a starting  %%%
%%% point for your own code with Stan, MATLABStan, and matstanlib.      %%%
%%%                                                                     %%%
%%% (c) beth baribault 2019 ---                            > matstanlib %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%% a bare bones outline of a Bayesian model fitting script
% 
% ...
% 
% to run this script, the following must be installed: 
%       Stan        >>  https://mc-stan.org/users/interfaces/cmdstan.html
%       MATLABStan  >>  https://github.com/brian-lau/MatlabStan
%       matstanlib  >>  https://github.com/baribault/matstanlib
%
% (this script was built from matstanlib's skeleton.m file.)

close all
clear
clc

%% inputs

modelName = 'CEST_2_pool';

%sampler settings
nChains     = 4;        %how many chains?
nWarmup     = 1000;     %how many iterations during warmup?
nIterations = 2000;     %how many iterations after warmup?

%an absolute path to a location for saved model output, figures, etc.:
% outputDir = 'X:\\my\custom\dir\'; %%% windows
% outputDir = '~/Documents/weizmann/rotation-1/CEST-fitting/';   %%% macOS, linux
outputDir = ['C:\Users\berna\Documents\Weizmann\rotation-1\CEST-fitting' filesep 'stan_output_' modelName];

%an absolute path to a (temporary) location for Stan output files:
workingDir = ['C:\Users\berna\Documents\Weizmann\rotation-1\CEST-fitting' filesep 'wdir_' modelName];

%delete temporary files after Stan finishes and output has been returned?
%(this helps prevent compilation conflicts, etc.)
cleanUp = true;

%% setup
fprintf('\n\n************\n\npreparing to run the %s model.\n',modelName)

%ensure all directories inputs end in a file separator
if ~isequal(outputDir(end),filesep), outputDir(end+1) = filesep; end
if ~isequal(workingDir(end),filesep), workingDir(end+1) = filesep; end

%output directory
fprintf(['\nmodel output and figures will be saved to the output directory:\n' ...
    '> %s\n'], outputDir)
if ~exist(outputDir,'dir')
    fprintf(['    currently this directory does not exist.  \n' ...
        '    making the directory now ... '])
    mkdir(outputDir)
    fprintf('done.\n')
end

%working directory
fprintf(['\nStan files will be saved to the working directory:\n' ...
    '> %s\n'], workingDir)
if ~exist(workingDir,'dir')
    fprintf(['    currently this directory does not exist.  \n' ...
        '    making the directory now ... '])
    mkdir(workingDir)
    fprintf('done.\n')
end

%% data
disp('simulating data ... ')
% disp('\nloading data ... ')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% either simulate or load data here %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clearvars gamma B0 w0 dendrite SEI electrolyte sys N

gamma = 16.546; % MHz/T; gyromagnetic ratio
B0 = 9.4; % T
w0 = gamma*B0; % in MHz
N = 200; % Number of measurements.

dendrite = struct('T1',1/8,'T2',1/393);
SEI = struct('T1',100,'T2',1/(28e3),'k',285,'f',0.02,'dw',-260*w0);
sys = struct('offsets',linspace(-500,500,N)*w0,'tp',0.2,'w1',500);

Z = CEST_multipool(sys,dendrite,SEI) + 0.01*randn(length(sys.offsets),1);
plot(Z,'.')

% Make sure that data is in [0,1]. Not sure why that is in the first place.
xZ = sys.offsets'; Z(Z >= 1) = 0.9995; Z(Z <= 0) = 1e-4;
% Save simulation data to csv format
writematrix([xZ,Z],'LP30_323.csv')

% collect data for Stan
dataStruct = struct('N',N,'xZ',sys.offsets,'Z',Z,'w1',sys.w1,'tp',sys.tp,'sigma',std(Z));

disp('done!')
%% Read actual data from table
gamma = 16.546; % MHz/T; gyromagnetic ratio of Li
B0 = 9.4; % T
w0 = gamma*B0; % in MHz

cd 'C:\Users\berna\Documents\Weizmann\rotation-1\CEST-fitting\data\derived';

T = readtable('LP30_dendrotes_CEST_exp_fit.xlsx',...
    'Range','A4:AA55');

% Make sure that data is in [0,1]. Not sure why that is in the first place.
xZ = T.Var1 * w0; Z = T.Var20;
Z(Z >= 1) = 0.9995;

% Save current experiment as csv, with name
% '@Material_@temperature_@RF_amplitude.csv'
writematrix([xZ,Z],'LP30_323_500.csv')

%% model specification

modelCode = {
    'functions { '
    '   // function declarations and definitions'
    '   real CEST(real offset, real T1a, real T2a, real T1b, real T2b, real f, real k, real dw, real t, real w1) {'
    '       real ka;    '
    '       vector[6] M;  '
    '       matrix[6,6] A;  '
    '       vector[6] b;  '
    ''
    '       b = [0, 0, 1/T1a, 0, 0, f/T1b]; '
    '       A = [ [-(1/T2a + ka), offset, 0, k, 0, 0]                   '
    '             [-offset, -(1/T2a + ka), w1, 0, k, 0]                 '
    '             [0, -w1, -(1/T1a + ka), 0, 0, k]                      '
    '             [ka, 0, 0, -(1/T2b + k), offset + dw*2*pi(), 0]       '
    '             [0, ka, 0, -(offset + dw*2*pi()), -(1/T2b + k), w1]   '
    '             [0, 0, ka, 0, -w1, -(1/T1b + k)] ];                   '
    '       M = scale_matrix_exp_multiply(t, A, b + A\b) - A\b;    '
    '       return M[3];'
    '   }'
    '}'
    'data { '
    '   int<lower=0> N; // number of measurements   '
    '   real xZ[N]; // offsets measured in spectrum   '
    '   real Z[N]; // experimental Z spectrum data    '
    '   real<lower=0> w1; // RF field intensity (in Hz)   '
    '   real<lower=0> tp; // saturation time (in s)'
    '   real sigma; // standard deviation of Z data '
    '}'
    ''
    'transformed data { '
    '  // transformations of data and fixed variables are declared here.'
    '  // convert offsets and RF field to rad/s'
    '   xZ = xZ*2*pi();   '
    '   w1 = w1*2*pi();   '
    '}'
    ''
    'parameters { '
    '   real<lower=0> T1a;  '
    '   real<lower=0> T2a;  '
    '   real<lower=0> T1b;  '
    '   real<lower=0> T2b;  '
    '   real<lower=0, upper=1> f;  '
    '   real<lower=0> k;  '
    '   real<upper=0> dw;    '
    '}'
    ''
    'transformed parameters {'
    '  //transformations and derived parameter declarations go here.'
    '  //...'
    '}'
    ''
    'model { '
    '  //priors and likelihood go here.'
    ''
    '   real Zmean[N];          '
    '   T1a ~ normal(8,0.8);    '
    '   T2a ~ normal(393,5);    '
    '   T1b ~ uniform(1,100);   '
    '   T2b ~ normal(28000,150); '
    '   f ~ uniform(0,1);   '
    '   k ~ normal(285,5);  '
    '   dw ~ normal(-260,2);  '
    ''
    '   for(i in 1:N){    '
    '       Zmean[i] = CEST(xZ[i],T1a,T2a,T1b,T2b,f,k,dw,tp,w1);'
    '       Z[i] ~ normal(Zmean[i],sigma);'
    '   }'
    '}'
    ''
    'generated quantities {'
    '  //predictive distributions and other quantites to track go here.'
    '  //...'
    '}'
};

%write the model code to a .stan file
stanFilePath = writestanfile(modelCode,modelName,workingDir);

%% compile the model
tic
fprintf('\ncompiling the model ... ')
sm = StanModel('file',stanFilePath);
sm.compile();
fprintf('done!\n')
fprintf('compiling took %.2f seconds.\n',toc)

%% run the model
tic
fprintf('\nrunning the model ... \n\n************\n\n')
fit =  sm.sampling('file',          stanFilePath, ...
                   'model_name',    modelName, ...
                   'sample_file',   modelName, ...
                   'verbose',       true, ...
                   'chains',        nChains, ...
                   'warmup',        nWarmup, ...
                   'iter',          nIterations, ...
                   'data',          dataStruct, ...
                   'working_dir',   workingDir);
fit.block();
[stanSummaryTxt,stanSummary] = fit.print('sig_figs',5);
fprintf('\n************\n\ndone!\n\n')

runtime = datevec(seconds(toc));
fprintf(['sampling took %i days, %i hours, %i minutes, ' ...
    'and %.2f seconds.\n'],runtime(3:end)) %optimistically assuming your 
                                           %model takes < 1 month to run :)

%% extract from the StanFit object
[samples,diagnostics] = extractsamples('MATLABStan',fit);
parameters = fieldnames(samples);
instances = getparaminstances([],samples);

%clean up after MATLABStan
clearvars fit
%clean up after Stan
if cleanUp, delete([workingDir '*']); rmdir(workingDir); end

%% save output
save([pwd filesep modelName '.mat'])

%% diagnostic reports & plots
%compute posterior sample-based diagnostics and summary statistics
posteriorTable = mcmctable(samples);
%print a report about all MCMC diagnostics
interpretdiagnostics(diagnostics,posteriorTable)

%trace plots/rank plots
tracedensity(samples,'???',diagnostics)
rankplots(samples,'???')

%% parameter estimates & model comparison & other statistics
estimatedValues = getsamplestats(samples);
% estimatedValues = getsamplestats(samples,trueValues);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%   call plotting functions, etc.   %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% other analyses

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%   compute other statistics here   %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% generate plots

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%   call plotting functions, etc.   %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
