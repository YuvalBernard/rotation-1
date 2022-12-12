function [Mz,Mx,My] = bloch(T1,T2,w1,dw,t0,tmax,q)
% Simulation of Bloch equations by solving differential equations in the
% rotating frame. w1 is defined is gamma*B1 and represents the power of the
% RF field. dw is the offset from resonance = gamma*B0 + w, where w is the
% RF field frequency.

% Assuming that B0 is initialy aligned on z-axis and B1 aligned on x-axis.
% In the rotating frame we have a system of ODEs:
% dMx/dt = dw*My(t) - Mx(t)/T2
% dMy/dt = -dw*Mx(t) - My(t)/T2 + w1*Mz(t)
% dMz/dt = -w1*My(t) - [Mz(t)-Mz0]/T1

% System of equations can be reduced to the form:
% M(t) = exp(-K*t)*[M0 - M_ss] + M_ss
% exp(-K*t) can be expressed as V*exp(-D*t)/V
% where [V,D] = eig(K).

% According to Cramer's method, M{i}_ss can be calculated by substituting
% the vector 'b' in the {i}th column of 'K' (representing {i} magnetization)
% and dividing between the determinants of the new and original 'K's.

Mz0 = 1; % initial z-magnetization
t = t0:(tmax-t0)/(q-1):tmax;
% Define propagation matrix
K = [1/T2 -dw 0;...
    dw 1/T2 -w1;...
    0 w1 1/T1];
% Define initial magnetization vector
b = [0 0 Mz0/T1]';
% Calculate solid-state values.
Mz_ss = det([K(:,1:2) b])/det(K);
Mx_ss = det([b K(:,2:3)])/det(K);
My_ss = det([K(:,1) b K(:,3)])/det(K);
M_ss = [Mx_ss My_ss Mz_ss]';

M = zeros(3,q);
for i = 1:q
    M(:,i) = fastExpm(-K*t(i))*(T1*b - M_ss) + M_ss;
end
Mz = M(3,:);

% Calculate also t_ss
if nargout >= 3
    t_ss = bloch_t_ss(t,Mz,Mz_ss);
    My = M(2,:);
    Mx = M(1,:);
end


