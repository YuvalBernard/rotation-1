function [Z,A,domain] = CEST_system_control(T1a,T2a,T1b,T2b,kb,M0a,M0b,dwa,db,w1,varargin)
% System control model of Bloch-McConnell equations to calculate Z and
% MTR_Assym Spectra
% Can also be used to find saturtion times of pool a and pool b. [Not
% implemented yet]

% Define constants
R1a = 1/T1a; R2a = 1/T2a;
R1b = 1/T1b; R2b = 1/T2b;
ka = kb*M0b/M0a; % Given from mass equality on both sides of exchange chemical equation
Z = zeros(size(dwa));
if ~isempty(varargin) % Obtain RF pulse duration if given
    tp = varargin{1};
end

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
    % Construct state-space two-pool model
    K = [-(R2a+ka) dwa(j) 0 kb 0 0;
        -dwa(j) -(R2a+ka) w1 0 kb 0;
        0 -w1 -(R1a+ka) 0 0 kb;
        ka 0 0 -(R2b+kb) (dwa(j)+db) 0;
        0 ka 0 -(dwa(j)+db) -(R2b+kb) w1;
        0 0 ka 0 -w1 -(R1b+kb)];

    B = [0 0 R1a*M0a/M0b 0 0 R1b]';
    C = [0 0 1 0 0 0];
    D = 0;
    sys = ss(K,B,C,D);
    if isempty(varargin) % Full saturation assumed
        T = M0b*tf(sys);
        Mza_ss = dcgain(T);
        Z(j) = Mza_ss/M0a;
    else
        t = 0:9e-4:tp;
        u = M0b*ones(size(t));
        x0 = [0 0 M0a 0 0 M0b]';
        Mza = lsim(sys,u,t,x0);
        Z(j) = Mza(end)/M0a;
    end
    
    if dwa(end) ~= -dwa(1) && j >= idx && j <= idx+domain
        K_ = -[-(R2a+ka) -dwa(j) 0 kb 0 0;
            dwa(j) -(R2a+ka) w1 0 kb 0;
            0 -w1 -(R1a+ka) 0 0 kb;
            ka 0 0 -(R2b+kb) (-dwa(j)+db) 0;
            0 ka 0 -(-dwa(j)+db) -(R2b+kb) w1;
            0 0 ka 0 -w1 -(R1b+kb)];

        sys = ss(K_,B,C,D);
        if isempty(varargin)
            T = M0b*tf(sys);
            Mza_ss_ = dcgain(T);
            A(j) = (Mza_ss_ - Mza_ss)/M0a;
        else
            t = 0:9e-4:tp;
            Mza_ = lsim(sys,u,t,x0);
            A(j) = (Mza_(end) - Mza(end))/Ma0;
        end
    end
end
if dwa(end) == -dwa(1)
    %  Assuming the domain is symmetric
    domain = round(length(dwa)/2);
    A = fliplr(Z(domain:end)) - Z(1:domain);
end
