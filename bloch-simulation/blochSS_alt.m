function [Mz_ss, t_ss] = blochSS_alt(varargin)

% This function gets the z-magnetization and calculates the magnetization
% in steady state and the time needed to reach it.
% Fit spline to data and calculate derivative.

T1 = varargin{1};
T2 = varargin{2};
w1 = varargin{3};
dw = varargin{4};

%%%%%%%%%%%%%%% Calculating Mz_ss %%%%%%%%%%%%%%%
Mz0 = 1; % initial z-magnetization
% Define propagation matrix
K = [1/T2 -dw 0;...
    dw 1/T2 -w1;...
    0 w1 1/T1];
% Define initial magnetization vector
b = [0 0 Mz0/T1]';
% Calculate solid-state value.
Mz_ss = det([K(:,1:2) b])/det(K);

%%%%%%%%%%%%%%% Calculating t_ss %%%%%%%%%%%%%%%
if nargout > 1 % calculate only if requested
    t0 = varargin{5}; tmax = varargin{6}; q = varargin{7};
    Mz = varargin{8};
    t = t0:(tmax-t0)/(q-1):tmax;
    % fit data to cubic spline
    f = spline(t,Mz);
    % evaluate derivatives at all t
    dMz = fnval(fnder(f,t));
    % Find t such that dMz -> 0, arbitrarily dMz < 1e-4.
    % Find the 50 first occurences that satisfy condition.
    poss_t_ss = find(abs(dMz) < 1e-4,50,'first');
    % If all occurences are adjacent, t_ss shall be the first occurence.
    % If the condition is not met [see T1,2=1,0.005; w1=150,dw=0 for reference]
    % take the first element of the subseries for which the occurences
    % are adjacent until the end
    j = 1;
    for i=1:(length(poss_t_ss)-1)
        if (poss_t_ss(i+1) == poss_t_ss(i) + 1)
            continue
        else
            j = i+1;
        end
    end
    t_ss = t(poss_t_ss(j));
end