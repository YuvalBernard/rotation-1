function [Mz,Mx,My] = bloch(t,T1,T2,w1,dw)
% Simulation of Bloch equations by solving differential equations in the
% rotating frame. w1 is defined is gamma*H1 and represents the power of the
% RF field. dw is the offset from resonance = gamma*H0 + w, where w is the
% RF field frequency.

% Define constants
R1 = 1/T1; % s^-1
R2 = 1/T2; % s^-1
Mz0 = 1; % initial z-magnetization

% Note: Assuming that H0 is initialy aligned on z-axis and H1 aligned on x-axis. 
% In the rotating frame we have a system of ODEs:
% dMx/dt = dw*My(t) - R2*Mx(t)
% dMy/dt = -dw*Mx(t) - R2*My(t) + w1*Mz(t)
% dMz/dt = -w1*My(t) - R1[Mz(t)-Mz0]
 
% These can be reduced to the matrix form:
% dM/dt = A*M
 
% The solution is of the form:
% M(t) = M(0)*exp(A*t)
% exp(A*t) can be expressed as:
% exp(A*t) = T*exp(D*t)*(1/T)
 
 A = [-R2 dw 0 0;...
     -dw -R2 w1 0;...
     0 -w1 -R1 R1*Mz0;...
     0 0 0 0];
 
 M0 = [0;0;1;1];
% We need the exponential variant of A, where calculation is done elementwise
M = zeros(length(M0),length(t));
for i = 1:length(t)
    M(:,i) = expm(A*t(i))*M0;
end
Mx = M(1,:); My = M(2,:); Mz = M(3,:);