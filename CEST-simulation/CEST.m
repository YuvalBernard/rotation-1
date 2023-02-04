function [Z,A,domain] = CEST(T1a,T2a,T1b,T2b,kb,M0a,M0b,dwa,db,w1,varargin)
% Simulation of Bloch-McConnell equations
% by solving differential equations in the
% rotating frame.

% Basically bloch equations plus exchange
% between pool a (water/bulk) and pool b (labile H/solute).

% Output is the Z and A values corresponding
% to Z = Mza(dwa)/M0a and A = [Mza(-dwa) - Mza(dwa)]/M0a
% which can be used to plot Z and MTR_Assymetry spectra.

% Calculation of Z and A spectra is available for both full and partial
% saturation of pool b. Saturation level is calculated if given RF (CW)
% pulse duration.

% dMxa/dt = -(R2a + ka)*Mxa(t) + dwa*Mya(t) + kb*Mxb(t)
% dMya/dt = -dwa*Mxa(t) - (R2a + ka)*Mya(t) - w1*Mza(t) + kb*Myb(t)
% dMza/dt = w1*Mya(t) - (R1a + ka)*Mza(t) + kb*Mzb(t) + R1a*Mza0
% dMxb/dt = ka*Mxa(t) - (R2b + kb)*Mxb(t) + dwb*Myb(t)
% dMyb/dt = ka*Mya(t) - dwb*Mxb(t) - (R2b + kb)*Myb(t) - w1*Mzb(t)
% dMzb/dt = ka*Mza(t) + w1*Myb(t) - (R1b + kb)*Mzb(t) + R1b*Mz0b


% dw{a,b} is the offset of a,b from resonance
% w1 represents the RF field "power": gamma*B1
% R{1,2}{a,b} is the reciprocal of T{1,2}{a,b}
% ka is the exchange rate from spins in pool a to pool b
% kb is the exchange rate from spins in pool b to pool a

% The equations above can be reduced to the matrix form:
% dM/dt = -K*M + b;

% The solution is of the form:
% M(t) = expm(-K*t)*(M0 - M_ss) + M_ss;     M_ss = K\b
% M = [Mxa(t) Mya(t) Mza(t) Mxb(t) Myb(t) Mzb(t)]'

% Define constants
R1a = 1/T1a; R2a = 1/T2a;
R1b = 1/T1b; R2b = 1/T2b;
ka = kb*M0b/M0a; % Given from mass equality on both sides of exchange chemical equation
b = [0; 0; R1a*M0a; 0; 0; R1b*M0b];
if ~isempty(varargin) % Obtain RF pulse duration if given
    tp = varargin{1};
end

w1 = w1*2*pi;
dwa = dwa*2*pi;
db = db*2*pi;

M0 = [0; 0; M0a; 0; 0; M0b];
Z = zeros(size(dwa));
% if dwa is not symmetric wrt 0, calculate A via the same method of Z
if dwa(end) ~= -dwa(1)
    % Find boundry for length of A
    if db < 0 % Peak should be at postive dwa
        if dwa(1) < 0 % Domain should contain all positive dwa slots
            idx = find(dwa >= 0,1);
            domain = length(dwa) - idx;
        else
            domain = find(dwa < 0,1) -1; idx = 0;
            if isempty(domain); domain = length(dwa); end
        end
    else % Peak should be at negative dwa
        if dwa(1) < 0 % Domain should contain all negative dwa slots
            domain = find(dwa > 0,1) -1; idx = 0;
            if isempty(domain); domain = length(dwa); end
        else
            idx = find(dwa <= 0,1);
            domain = length(dwa) - idx;
        end
    end

    A = zeros(1,domain); % Initialize A
end

for j = 1:length(dwa)

    K = [(R2a+ka) -dwa(j) 0 -kb 0 0;
        dwa(j) (R2a+ka) -w1 0 -kb 0;
        0 w1 (R1a+ka) 0 0 -kb;
        -ka 0 0 (R2b+kb) -(dwa(j)+db) 0;
        0 -ka 0 (dwa(j)+db) (R2b+kb) -w1;
        0 0 -ka 0 w1 (R1b+kb)];

    M_ss = K\b;

    if isempty(varargin) % Full saturation assumed
        Z(j) = M_ss(3)/M0a;
    else % Calculate Z for given tp
        M_tp = fastExpm(-K*tp)*(M0 - M_ss) + M_ss;
        Z(j) = M_tp(3)/M0a;
    end

    if dwa(end) ~= -dwa(1) && j >= idx && j <= idx+domain

        K_ = [(R2a+ka) dwa(j) 0 -kb 0 0;
            -dwa(j) (R2a+ka) -w1 0 -kb 0;
            0 w1 (R1a+ka) 0 0 -kb;
            -ka 0 0 (R2b+kb) -(-dwa(j)+db) 0;
            0 -ka 0 (-dwa(j)+db) (R2b+kb) -w1;
            0 0 -ka 0 w1 (R1b+kb)];

        M_ss_ = K_\b;

        if isempty(varargin) % Full saturation assumed
            A(j) = (M_ss_(3) - M_ss(3))/M0a;
        else % Calculate A for given tp
            M_tp_ = fastExpm(-K_*tp)*(M0 - M_ss_) + M_ss_;
            A(j) = (M_tp_(3) - M_tp(3))/M0a;
        end
    end
end
if dwa(end) == -dwa(1)
    %  Assuming the domain is symmetric
    domain = round(length(dwa)/2);
    A = fliplr(Z(domain:end)) - Z(1:domain);
end
return