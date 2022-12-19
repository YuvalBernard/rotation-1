function [Z,A,domain] = CEST_alt(T1a,T2a,T1b,T2b,kb,M0a,M0b,dwa,db,w1,tp)

% Define constants
R1a = 1/T1a; R2a = 1/T2a;
R1b = 1/T1b; R2b = 1/T2b;
ka = kb*M0b/M0a; % Given from mass equality on both sides of exchange chemical equation
M0 = [0 0 0 0 M0a M0b 1]';
dwb = dwa + db;
Z = zeros(size(dwa));
for j = 1:length(dwa)
    A = [-(R2a+ka) kb dwa(j) 0 0 0 0;
        ka -(R2b+kb) 0 dwb(j) 0 0 0;
        -dwa(j) 0 -(R2a+ka) kb w1 0 0;
        0 -dwb(j) ka -(R2b+kb) 0 w1 0;
        0 0 -w1 0 -(R1a+ka) kb R1a*M0a;
        0 0 0 -w1 ka -(R1b+kb) R1b*M0b;
        0 0 0 0 0 0 0];
    Mprime = @(y) A*y;
    [t,y] = ode45(@(t,y) Mprime(y),[0 tp],M0);
%     M_tp = expm(A*tp)*M0;
%     Z(j) = M_tp(5)/M0a;
Z(j) = y(end,5)/M0a;
end
%  Assuming the domain is symmetric
domain = round(length(dwa)/2);
A = fliplr(Z(domain:end)) - Z(1:domain);
