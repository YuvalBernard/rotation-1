function [Mz, Mx, My, t] = bloch_numerical(T1,T2,w1,dw,t0,tmax,q)

A = [-1/T2 dw 0 0;...
     -dw -1/T2 w1 0;...
     0 -w1 -1/T1 1/T1;... 
     0 0 0 0];

Mprime_alt = @(y) A*y;
[t,y] = ode45(@(t,y) Mprime_alt(y),[t0 tmax],[0;0;1;1]);
% y = deval(sol,t0:(tmax-t0)/(q-1):tmax);

% Mx = y(1,:); My = y(2,:); Mz = y(3,:);
Mx = y(:,1); My = y(:,2); Mz = y(:,3);
