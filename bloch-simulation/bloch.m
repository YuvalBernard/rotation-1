function [Mz,Mx,My] = bloch(T1,T2,w1,dw,t0,tmax,q)
% Simulation of Bloch equations by solving differential equations in the
% rotating frame. w1 is defined is gamma*H1 and represents the power of the
% RF field. dw is the offset from resonance = gamma*H0 + w, where w is the
% RF field frequency.

% Assuming that H0 is initialy aligned on z-axis and H1 aligned on x-axis. 
% In the rotating frame we have a system of ODEs:
% dMx/dt = dw*My(t) - Mx(t)/T2
% dMy/dt = -dw*Mx(t) - My(t)/T2 + w1*Mz(t)
% dMz/dt = -w1*My(t) - [Mz(t)-Mz0]/T1
 
% These can be reduced to the matrix form:
% dM/dt = A*M
 
% The solution is of the form:
% M(t) = M(0)*exp(A*t)
% exp(A*t) can be expressed as V*exp(D*t)*(1/V)
% where [V,D] = eig(A).

Mz0 = 1; % initial z-magnetization

 A = [-1/T2 dw 0 0;...
     -dw -1/T2 w1 0;...
     0 -w1 -1/T1 Mz0/T1;...
     0 0 0 0];
 
 M0 = [0;0;1;1];
 
% Calculate M(t) by calculating expA.
% We need the exponential variant of A, where calculation is done elementwise
% Initialize M(t)
M = zeros(4,q);
% Calculate eigenvalues of A
% [V,D] = eig(A);
% Save computation time by preprocessing exp(diag(D))
% expdiag = exp(diag(D));
t = t0:(tmax-t0)/(q-1):tmax;
for i=1:q
%     M(:,i) = V*diag(expdiag.^t(i))/V*M0;
    M(:,i) = expm(A*t(i))*M0;
end
M = real(M);
Mx = M(1,:); My = M(2,:); Mz = M(3,:);
