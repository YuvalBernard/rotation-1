function [Z,MTR] = CEST_multipool(System,varargin)
% Attempt at solving multi-pool CEST of given dimensionality.
% Calculate Z spectrum and corresponding MTR Asymmetry spectrum
% Note: Equilibrium magnetization of large pool normalized to 1.
% Assume that each pool exchanges only with pool a!
% System: Struct containing:
% w1 1x1 (Hz) - RF amplitude
% offsets Nx1 (Hz) - measured range of offsets from resonance of pool a
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

offsets = System.offsets; N = length(offsets);
n = nargin-1; % Number of pools
Ex = zeros(n,n); % Coefficient matrix of system of linear equations representing exchange
S = zeros(3*n,3*n); % Block diagonal matrix containing the sufficient properties of each pool
% in order to solve bloch equations for it (T1,T2,dw,w1)
b = zeros(3*n,1); % vector from the system dM/dt = -A*M + b
M0 = zeros(size(b)); % vector containing initial magnetizations;

%%%%%%%%%%%%% Create necessary matrices to solve Bloch-McConnell equations
j = 1; % Auxiliary index
for i=1:n
    P_i = varargin{i}; % Pool i struct
    if i == 1; P_i.f = 1; P_i.dw = 0; end
    % Verify that pool a was given first.
    if i == 1 && isfield(P_i,'k'); error('Large pool not assigned first.'); end
    % Populate S
    if ~isfield(P_i,'T1'); error('T1 of pool %s not provided.', num2str(i)); end
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
end
K = S + kron(Ex,eye(3)); % Coefficient matrix of Bloch-McConnell equations of all pools.
% Only parameter missing is the saturation offset, which changes each interation.

%%%%%%%%%%%%% Calculate Z spectrum
% Now that the system parameters are almost fully initialized, calculate Z
% spectum.
Z = zeros(N,1);

for j = 1:N
    % Finish initialization of full coefficient matrix.
    A = K + kron(eye(n),[0 offsets(j) 0;
                        -offsets(j) 0 0;
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
end

%%%%%%%%%%%%% Calculate MTR spectrum
% MTR spectrum is calculated by taking Z_ref - Z,
% where Z_ref is the associated Z spectrum without exchange.
% Just solve BM equations again but with P_i.f set to zero.
% Reconstruct core matrices:
Ex(1:end,1) = 0; b(6:end) = 0; M0(6:end) = 0;
K = S + kron(Ex,eye(3));
% Initizlize MTR spectrum array
MTR = zeros(N,1);
for j = 1:N
    % Finish initialization of full coefficient matrix.
    A = K + kron(eye(n),[0 offsets(j) 0;
                        -offsets(j) 0 0;
                        0 0 0]);
    
    % Calculate steady-state magnetization according to Bloch McConnell
    % equations
    M_ss = A\b;
    if ~isfield(System,'tp') % Full saturation assumed
        MTR(j) = M_ss(3) - Z(j);
    else % Calculate Z for given tp
        M_tp = fastExpm(-A * System.tp)*(M0 - M_ss) + M_ss;
        MTR(j) = M_tp(3) - Z(j);
    end
end