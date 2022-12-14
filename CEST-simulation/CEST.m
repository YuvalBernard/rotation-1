function [Z,A,domain,t_ss] = CEST(T1a,T2a,T1b,T2b,kb,M0a,M0b,dwa,dwb,w1,varargin)
% Simulation of Bloch-McConnell equations
% by solving differential equations in the
% rotating frame.

% Basically bloch equations plus exchange
% between pool a (water/bulk) and pool b (labile H/solute).

% Output is the Z and A values corresponding
% to Z = Mza(dwa)/M0a and A = [Mzb(-dwa) - Mzb(dwa)]/M0a
% which can be used to plot Z and Assymetry spectra.

% Calculation of Z and A spectra is available for both full and partial
% saturation of pool b. Saturation level is calculated if given RF (CW)
% pulse duration.

% dMxa/dt = dwa*Mya(t) - R2a*Mxa(t) - ka*Mxa(t) + kb*Mxb(t)
% dMxb/dt = dwb*Myb(t) - R2b*Mxb(t) - kb*Mxb(t) + ka*Mxa(t)
% dMya/dt = -dwa*Mxa(t) - R2a*Mya(t) - ka*Mya(t) + kb*Myb(t) + w1*Mza(t)
% dMyb/dt = -dwb*Mxb(t) - R2b*Myb(t) - kb*Myb(t) + ka*Mya(t) + w1*Mzb(t)
% dMza/dt = -w1*Mya(t) - R1a*[Mza(t)-M0a] - ka*Mza(t) + kb*Mzb(t)
% dMzb/dt = -w1*Myb(t) - R1b*[Mzb(t)-M0b] - kb*Mzb(t) + ka*Mza(t)

% From paper:
% dMxa/dt = -(R2a + ka)*Mxa(t) - dwa*Mya(t) + kb*Mxb(t)
% dMya/dt = dwa*Mxa(t) - (R2a + ka)*Mya(t) - w1*Mza(t) + kb*Myb(t)
% dMza/dt = w1*Mya(t) - (R1a + ka)*Mza(t) + kb*Mzb(t) + R1a*Mza0
% dMxb/dt = ka*Mxa(t) - (R2b + kb)*Mxb(t) - dwb*Myb(t)
% dMyb/dt = ka*Mya(t) + dwb*Mxb(t) - (R2b + kb)*Myb(t) - w1*Mzb(t)
% dMzb/dt = ka*Mza(t) + w1*Myb(t) - (R1b + kb)*Mzb(t) + R1b*Mz0b


% dw{a,b} is the offset of a,b from resonance
% w1 represents the RF field "power": gamma*B1
% R{1,2}{a,b} is the reciprocal of T{1,2}{a,b}
% ka is the exchange rate from spins in pool a to pool b
% kb is the exchange rate from spins in pool b to pool a

% The equations above can be reduced to the matrix form:
% dM/dt = F*M
% or dM/dt = -K*M + b;

% The solution is of the form:
% M(t) = expm(F*t)*M0;
% M = [Mxa(t) Mya(t) Mza(t) Mxb(t) Myb(t) Mzb(t) 1]'
% or M(t) = expm(-K*t)*(M0 - M_ss) + M_ss;
% M = [Mxa(t) Mya(t) Mza(t) Mxb(t) Myb(t) Mzb(t)]'

% Define constants
R1a = 1/T1a; R2a = 1/T2a;
R1b = 1/T1b; R2b = 1/T2b;
ka = kb*M0b/M0a; % Given from mass equality on both sides of exchange chemical equation
if ~isempty(varargin) % Obtain RF pulse duration if given
    tp = varargin{1};
    t = t0:1e-6:tp;
end
M0 = [0 0 M0a 0 0 M0b]';
Z = zeros(size(dwa));


if nargout < 4 % No dynamics analysis (t_ss) requested.
    for j = 1:length(dwa)
       
        K = -[-(R2a+ka) -dwa(j) 0 kb 0 0;...
            dwa(j) -(R2a+ka) -w1 0 kb 0;...
            0 w1 -(R1a+ka) 0 0 kb;...
            ka 0 0 -(R2b+kb) -dwb(j) 0;...
            0 ka 0 dwb(j) -(R2b+kb) -w1;...
            0 0 ka 0 w1 -(R1b+kb)];

        b = [0 0 R1a*M0a 0 0 R1b*M0b]';
        M_ss = K\b;
        if isempty(varargin) % Full saturation assumed
            Z(j) = M_ss(3)/M0a;
        else % Partial saturation possible
            M_tp = fastExpm(-K*tp)*(M0 - M_ss) + M_ss;
            Z(j) = M_tp(3)/M0a;
        end
    end
    domain = round(length(dwa)/2);
    A = fliplr(Z(domain:end)) - Z(1:domain);
    return
end
% We're asked to calculate t_ss. return table of t_sat and t_ss
t_ssa = zeros(size(dwb));
tsatb = zeros(size(dwb));

t_ss = table(t_satb,t_ssa,'VariableNames',{'t^a_{ss}','t^b_{sat}'});

