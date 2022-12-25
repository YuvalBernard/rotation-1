function t_ss = calc_t_ss(T1,T2,w1,dw)
% Different approaches to calculate time required to approach steady state

%%%%%%%%%%% Disclaimer: %%%%%%%%%%%
% I don't think it can be found analytically.
% Also, Approach 1 seems to be more consistent than Approach 2.

%%%%%%%%%%% Approach 1  %%%%%%%%%%%
% Stability analysis dictates that:
% Given a system of linear ODEs such that dx/dt = A*x + b
% (a) There is a unique solution if and only if rank(A) = rank(A|b)
% (b) Solution is stable if the real parts of the eigenvalues of A are negative
% (c) The Solution decays exponentially, relatively to the maximal (negative) eigenvalue of A.

% Thus, the settling time can be approximated to desired accuracy as follows:
% t = x/EV
% x = ln(100/error), error in percentage (represents desired accuracy)
% EV = min(abs(real(eig(A))))

A = [1/T2 -dw 0;...
    dw 1/T2 -w1;...
    0 w1 1/T1];

err = 0.005; % Set this to desired degree of accuracy
t_ss = log(100/err)/real(eigs(A,1,'smallestreal'));

%%%%%%%%%%% Approach 2  %%%%%%%%%%%
% Contstruct a state-space model of the system (LTI system)
% Extract t_ss using given tools provided by MATLAB's Control System
% Toolbox

% Differential equation for state variables:
% dx/dt = A*x + B*u
% where x = [Mx My Mz]', A is the coefficient matrix of the linear ODE, and
% B is [0 0 1/T1]'. u is the single input: Mz0.
% Equation for output:
% y = C*x + D*u
% where y = Mz is the output, C = [0 0 1], D = 0.

% A = [-1/T2 dw 0;
%     -dw -1/T2 w1;
%     0 -w1 -1/T1];
% B = [0 0 1/T1]';
% C = [0 0 1];
% D = 0;
% sys = ss(A,B,C,D);
% ST = 5e-4; % (%error/100)
% S = stepinfo(sys,'SettlingTimeThreshold',ST);
% t_ss = S.SettlingTime;

% SettlingTime is taken as the first time T such that
% error|y(t)-y_ss| <= ST*|y_ss-y(0)|