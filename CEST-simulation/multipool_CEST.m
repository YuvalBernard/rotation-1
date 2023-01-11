function [Z,Asym,domain] = multipool_CEST(System,varargin)
% Attempt at solving multi-pool CEST of given dimensionality.
% Note: Asymmetry spectra should be implemented if this works.
% Note: Here equilibrium z magnetization of pool a is normalized to 1.
% Assume that each pool exchanges only with pool a!
% System: Struct containing:
% w1 1x1 (Hz) - RF amplitude
% dwa Nx1 (Hz) - measured range of offsets from resonance of pool a
% tp 1x1 (s) - saturation time (optional)

% Pool_a is the pool of the free agent (e.g. water hydrogens)
% Pool_i is the pool of the i_th agent (e.g. metabolyte proton 1,2,3 etc.)
% pool i=1 is pool a.
% Each pool is a struct array containing:
% T1 1x1 (s) - longitudinal relaxation
% T2 1x1 (s) - transverse relaxation
% dw 1x1 (Hz) - offset from pool a
% f 1x1 (arb) - (effective) molar fraction compared to pool a
% k 1x1 (Hz) - chemical exchange rate with pool a

dwa = System.dwa;
n = nargin-1; % Number of pools
Ex = zeros(n,n); % Coefficient matrix of system of linear equations representing exchange
S = zeros(3*n,3*n); % Block diagonal matrix containing the sufficient properties of each pool
% in order to solve bloch equations for it (T1,T2,dw,w1)
b = zeros(3*n,1); % vector from the system dM/dt = -K*M + b
M0 = zeros(size(b)); % vector containing initial magnetizations;


%%%%%%%%%%%%% Create necessary matrices to solve Bloch-McConnell equations
j = 1; % Auxiliary index
if dwa(end) ~= -dwa(1); flag = 0; end % Flag that tracks whether peaks in Z spectrum are located on positive dwa, negative, or both.
for i=1:n
    P_i = varargin{i}; % Pool i struct
    if i == 1; P_i.f = 1; P_i.dw = 0; end
    % Verify that pool a was given first.
    if i == 1 && isfield(P_i,'k'); error('Large pool not assigned first.'); end
    % Populate S
    if ~isfield(P_i,'T1'); error('nT1 of pool %s not provided.', num2str(i)); end
    if ~isfield(P_i,'T2'); error('T2 of pool %s not provided.', num2str(i)); end
    if ~isfield(P_i,'dw'); error('Delta omega (offset) of pool %s not provided.', num2str(i)); end
    if ~isfield(System,'w1'); error('w1 (RF field amplitude) not provided.'); end
    S(j,j) = 1/P_i.T2;
    S(j,j+1) = P_i.dw;
    S(j+1,j) = -S(j,j+1);
    S(j+1,j+1) = S(j,j);
    S(j+1,j+2) = System.w1;
    S(j+2,j+1) = -System.w1;
    S(j+2,j+2) = 1/P_i.T1;
    j = j+3;
    % Populate b and M0
    b(j-1) = (1/P_i.T1) * P_i.f;
    M0(j-1) = P_i.f;
    % Populate Ex
    if i ==1; continue; end
    if ~isfield(P_i,'k'); error('Exchange rate with pool %s not provided.', num2str(i)); end
    if ~isfield(P_i,'f'); error('Fraction of pool %s not provided.', num2str(i)); end
    Ex(1,1) = Ex(1,1) + P_i.f * P_i.k;
    Ex(1,i) = -P_i.k;
    Ex(i,1) = -P_i.f * P_i.k;
    Ex(i,i) = P_i.k;

    if dwa(end) ~= -dwa(1)
        if P_i.dw < 0 % peak shoud be at positive dwa
            flag = 1;
        else % peak should be at negative dwa
            flag = - 1;
        end
    end
end
K = S + kron(Ex,eye(3)); % Coefficient matrix of Bloch-McConnell equations of all pools.
% Only parameter missing is dwa, which changes each interation.

%%%%%%%%%%%%% Select domain where Asym spectrum is defined.
% If dwa is not symmetric around zero, find domain of definition
% We assume here that both exchangable spins are located on the same side of the spectrum
% (both are either at negative or at positive ppms). Otherwise the Asym spectrum doesn't make sense...

if dwa(end) ~= -dwa(1)
    % Find boundry for length of Asym
    if flag > 0 % Peaks are at postive dwa
        if dwa(1) < 0 % Domain should contain all positive dwa slots
            idx = find(dwa >= 0,1);
            domain = length(dwa) - idx;
        else
            domain = find(dwa < 0,1) -1; idx = 0;
            if isempty(domain); domain = length(dwa); end
        end
    else % Peaks are at negative dwa
        if dwa(1) < 0 % Domain should contain all negative dwa slots
            domain = find(dwa > 0,1) -1; idx = 0;
            if isempty(domain); domain = length(dwa); end
        else
            idx = find(dwa <= 0,1);
            domain = length(dwa) - idx;
        end
    end

    Asym = zeros(domain,1); % Initialize Asym
end

%%%%%%%%%%%%% Calculate Z and Asym spectra
% Now that the system parameters are almost fully initialized, calculate Z
% spectum.
Z = zeros(size(dwa));

for j = 1:length(dwa)
    % Finish initialization of full coefficient matrix.
    A = K + kron(eye(n),[0 dwa(j) 0;
        -dwa(j) 0 0;
        0 0 0]);

    % Calculate steady-state magnetization according to Bloch McConnell
    % equations
    M_ss = A\b;
    if ~isfield(System,'tp') % Full saturation assumed
        Z(j) = M_ss(3);
    else % Calculate Z for given tp
        M_tp = fastExpm(-A * System.tp)*(M0 - M_ss) + M_ss;
        Z(j) = M_tp(3);
    end

    if dwa(end) ~= -dwa(1) && j >= idx && j <= idx+domain
        A = K - kron(eye(n),[0 dwa(j) 0;
            -dwa(j) 0 0;
            0 0 0]);

        M_ss_ = A_\b;

        if ~isfield(System,'tp') % Full saturation assumed
            Asym(j) = M_ss_(3) - M_ss(3);
        else
            M_tp_ = fastExpm(-A * System.tp)*(M0 - M_ss_) + M_ss_;
            Asym(j) = M_tp_(3) - M_tp(3);
        end
    end
    if dwa(end) == -dwa(1) % Domain must be symmetric
        domain = round(length(dwa)/2);
        Asym = fliplr(Z(domain:end)) - Z(1:domain);
    end
end
