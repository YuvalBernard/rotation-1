function [t_ss,Mz_ss] = blochSS_alt(T1,T2,w1,dw,t0,tmax,q,Mz)
% Different approach to calculating z-magnetization at steady state
% and time required to reach steady state.

% Mz_ss is calculated analytically using Cramer's method.
% Define propagation matrix:
t = t0:(tmax-t0)/(q-1):tmax;
Mz0 = 1; % initial z-magnetization

K = [1/T2 -dw 0;...
     dw 1/T2 -w1;...
     0 w1 1/T1];
% Define initial magnetization vector
b = [0 0 Mz0/T1]';

% bloch equations state that dM/dt = -K*M + b, M = [Mx My Mz]
% In steady state dM/dt = 0. In other words, K*M = b
% According to Cramer's method, Mz_ss can be calculated by substituting
% the vector 'b' in the third column of 'K' (representing z magnetization)
% And dividing between the determinants of the new and original 'K's.
Mz_ss = det([K(:,1:2) b])/det(K);

% Now we need to find t_ss. I don't think it can be found analytically
% Also, extrapolation is strongly unadvised.
%%% Approach 1: fit Mz-Mz_ss to spline and use fnzeros to find roots of 
%%% fit. Then, take the min root found.
% 
% f = spline(t,Mz-Mz_ss);
% t_ss = min(fnzeros(f),[],'all');

%%% Approach 2: find when Mz = Mz_ss and dMz/dt = 0 (minimal value).

idx = find((Mz <= Mz_ss + 5e-5) & (abs(gradient(Mz,t)) >= 1e-5 ),1);
t_ss = t(idx);

%%% Approach 3: numerically solve the following equation derived from the
%%% formal solution to bloch equations.
%%% M(t) = exp(-K*t)*(M0 - M_ss) + M_ss.
%%% Take only the part that applies to Mz at SS (and subtract Mz_ss from both sides):
%%% 0 = -Mx*ss*exp(-K(3,1)*t) - My_ss*exp(-K(3,2)*t) + (Mz0/T1-Mz_ss)*exp(-K(3,3)*t)

% Mx_ss = det([b K(:,2:end)])/det(K);
% My_ss = det([K(:,1) b K(:,3)])/det(K);
% syms x; % representing t_ss
% eqn = -Mx_ss*exp(-K(3,1)*x) - My_ss*exp(-K(3,2)*x) + exp(-K(3,3)*x)*(Mz0/T1 - Mz_ss) == 0;
% t_ss = double(vpasolve(eqn));