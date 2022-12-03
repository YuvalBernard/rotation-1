function [Mz, Mx, My] = bloch_numerical_alt(T1,T2,w1,dw,t0,tmax,q)

Mprime = @(y) [dw*y(2) - y(1)/T2;...
              -dw*y(1) - y(2)/T2 + w1*y(3);...
              -w1*y(2) - y(3)/T1 + 1/T1];

sol = ode45(@(t,y) Mprime(y),[t0 tmax],[0;0;1]);

y = deval(sol,t0:(tmax-t0)/(q-1):tmax);

Mx = y(1,:); My = y(2,:); Mz = y(3,:);