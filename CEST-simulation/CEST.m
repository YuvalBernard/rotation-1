function [Z,A,domain,t_ss] = CEST(T1a,T2a,T1b,T2b,kb,M0a,M0b,dwa,dwb,w1,t0,tmax,q)
% Simulation of Bloch-McConnell equations
% by solving differential equations in the
% rotating frame.

% Basically bloch equations plus exchange
% between pool a (water/bulk) and pool b (labile H/solute).

% Output is the Z and A values corresponding
% to Z = Mza(dwa)/M0a and A = [Mzb(-dwa) - Mzb(dwa)]/M0a
% which can be used to plot Z and Assymetry spectra.

% dMxa/dt = dwa*Mya(t) - R2a*Mxa(t) - ka*Mxa(t) + kb*Mxb(t)
% dMxb/dt = dwb*Myb(t) - R2b*Mxb(t) - kb*Mxb(t) + ka*Mxa(t)
% dMya/dt = -dwa*Mxa(t) - R2a*Mya(t) - ka*Mya(t) + kb*Myb(t) + w1*Mza(t)
% dMyb/dt = -dwb*Mxb(t) - R2b*Myb(t) - kb*Myb(t) + ka*Mya(t) + w1*Mzb(t)
% dMza/dt = -w1*Mya(t) - R1a*[Mza(t)-M0a] - ka*Mza(t) + kb*Mzb(t)
% dMzb/dt = -w1*Myb(t) - R1b*[Mzb(t)-M0b] - kb*Mzb(t) + ka*Mza(t)

% dw{a,b} is the offset of a,b from resonance
% w1 represents the RF field "power": gamma*B1
% R{1,2}{a,b} is the reciprocal of T{1,2}{a,b}
% ka is the exchange rate from spins in pool a to pool b
% kb is the exchange rate from spins in pool b to pool a

% The equations above can be reduced to the matrix form:
% dM/dt = F*M
 
% The solution is of the form:
% M(t) = M(0)*exp(F*t)
% M = [Mxa(t) Mxb(t) Mya(t) Myb(t) Mza(t) Mzb(t) 1]'

R1a = 1/T1a; R2a = 1/T2a;
R1b = 1/T1b; R2b = 1/T2b;
ka = kb*M0b/M0a;

M0 = [0; 0; 0; 0; M0a; M0b; 1];
M = zeros(length(M0),q);
t = t0:(tmax-t0)/(q-1):tmax;
Z = zeros(size(dwa));

if nargout < 4 % No dynamics requested. Calculate Z and A directly.
    for j = 1:length(dwa)
        K = [(R2a+ka) -kb -dwa(j) 0 0 0;...
            -ka (R2b+kb) 0 -dwb(j) 0 0;...
            dwa(j) 0 -(-R2a+ka) -kb -w1 0;...
            0 dwb(j) -ka (R2b+kb) 0 -w1;...
            0 0 w1 0 (R1a+ka) -kb;...
            0 0 0 w1 -ka (R1b+kb)];

        b = [0 0 0 0 R1a*M0a R1b*M0b]';
        Z(j) = (det([K(:,1:4) b K(:,6)])/det(K))/M0a;
    end
    domain = round(length(dwa)/2);
    A = fliplr(Z(domain:end)) - Z(1:domain);
    return
end

for j = 1:length(dwb)

    F = [-(R2a+ka) kb dwa(j) 0 0 0 0;...
        ka -(R2b+kb) 0 dwb(j) 0 0 0;...
        -dwa(j) 0 (-R2a+ka) kb w1 0 0;...
        0 -dwb(j) ka -(R2b+kb) 0 w1 0;...
        0 0 -w1 0 -(R1a+ka) kb R1a*M0a;...
        0 0 0 -w1 ka -(R1b+kb) R1b*M0b;...
        0 0 0 0 0 0 0];

    for k=1:q
%             M(:,k) = expm(F*t(k))*M0;
        M(:,k) = fastExpm(F*t(k))*M0;
    end
    f = spline(t,M(5,:)); % model Mza as 'spline'
    dfdt = fnder(f);
    dMzdt = fnval(dfdt,t);
    poss_t_ss = find(abs(dMzdt) < 1e-4,100,'first');
    jj = 1;
    for i=1:(length(poss_t_ss)-1)
        if (poss_t_ss(i+1) == poss_t_ss(i) + 1)
            continue
        else
            jj = i+1;
        end
    end
    t_ss = t(poss_t_ss(jj));
    Mza = fnval(f,t_ss);
    Z(j) = Mza/M0a;
end

domain = round(length(dwa)/2);
A = fliplr(Z(domain:end)) - Z(1:domain);

