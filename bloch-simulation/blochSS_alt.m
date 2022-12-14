function [Mz_ss, t_ss] = blochSS_alt(T1,T2,w1,dw)

% This function gets the z-magnetization and calculates the magnetization
% in steady state and the time needed to reach it.
% Fit spline to data and calculate derivative.

%%%%%%%%%%%%%%% Calculating Mz_ss %%%%%%%%%%%%%%%
Mz0 = 1; % initial z-magnetization
% Define propagation matrix
K = [1/T2 -dw 0;...
    dw 1/T2 -w1;...
    0 w1 1/T1];
% Define initial magnetization vector
b = [0 0 Mz0/T1]';
% Calculate solid-state value.
M_ss = K\b; Mz_ss = M_ss(3);

%%%%%%%%%%%%%%% Calculating t_ss %%%%%%%%%%%%%%%
if nargout > 1 % calculate only if requested
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

    A = [-1/T2 dw 0;...
        -dw -1/T2 w1;...
        0 -w1 -1/T1];

    EV = min(abs(real(eig(A))));
    err = 0.067; % Set this to desired degree of accuracy
    t_ss = log(100/err)/EV;

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
end