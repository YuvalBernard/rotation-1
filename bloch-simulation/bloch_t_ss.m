function t_ss = bloch_t_ss(varargin)
% Different approach to calculating time required to approach steady state
t = varargin{1};
Mz = varargin{2};
Mz_ss = varargin{3};

%%%%%%%%%%% Disclaimer: %%%%%%%%%%%
% I don't think it can be found analytically.
% Also, extrapolation is strongly unadvised.

%%% original approach
f = spline(t,Mz);
% evaluate derivatives at all t
dMz = fnval(fnder(f),t);
poss_t_ss = find(abs(dMz) < 1e-4,100,'first');
%%% If all occurences are adjacent, t_ss shall be the first occurence.
%%% If the condition is not met [see T1,2=1,0.005; w1=150,dw=0 for reference]
%%% take the first element of the subseries for which the occurences
%%% are adjacent until the end
j = 1;
for i=1:(length(poss_t_ss)-1)
    if (poss_t_ss(i+1) == poss_t_ss(i) + 1)
        continue
    else
        j = i+1;
    end
end
t_ss = t(poss_t_ss(j));

%%% Approach 2: fit Mz-Mz_ss to spline and use fnzeros to find roots of 
%%% fit. Then, take the min root found.

% f = spline(t,Mz-Mz_ss);
% t_ss = min(fnzeros(f),[],'all');

%%% Approach 3: find when Mz = Mz_ss and dMz/dt = 0 (minimal value).
% dMz = fnval(fnder(spline(t,Mz)),t);
% if Mz_ss <= 1e-2
%     idx = find((Mz <= Mz_ss + 5e-5) & (abs(dMz) <= 5e-5),1);
% end
% idx = find((Mz <= Mz_ss + 5e-3) & (abs(dMz) <= 5e-4),1);
% t_ss = t(idx);

% poss_t_ss = find(abs(dMz) < 1e-4,100,'first');

%%% Approach 4: numerically solve the following equation derived from the
%%% formal solution to bloch equations.
%%% M(t) = exp(-K*t)*(M0 - M_ss) + M_ss.
%%% Take only the part that applies to Mz at SS (and subtract Mz_ss from both sides):
%%% 0 = -Mx*ss*exp(-K(3,1)*t) - My_ss*exp(-K(3,2)*t) + (Mz0/T1-Mz_ss)*exp(-K(3,3)*t)

% if nargin > 3
%     T1 = varargin{4};
%     M_ss = varargin{5};
%     % Define propagation matrix and initial z magnetization:
%     K = [1/T2 -dw 0;...
%         dw 1/T2 -w1;...
%         0 w1 1/T1];
%     Mz0 = 1;
% 
%     problem.objective = @(x) -M_ss(1)*exp(-K(3,1)*x) - M_ss(2)*exp(-K(3,2)*x) + exp(-K(3,3)*x)*(Mz0/T1 - M_ss(3));
%     problem.x0 = 0;
%     problem.solver = 'fzero';
%     oldopts = optimset('Display','iter','TolX',1e-4,'TolFun',5e-5); % set tolerance to 5e-5.
%     newopts = optimset('Display','iter');
%     problem.options = optimset(oldopts,newopts);
% 
%     t_ss = fzero(problem);
% end

%%% Approach 4: fit the inverse spline and find the value of spline(Mz,t) at Mz = Mz_ss.
%%% First, average data points with the same site
% [x,y] = chckxywp(Mz,t);
% f = spline(x,y);
% idx = find(x >= Mz_ss + 1e-4,1,"first");
% t_ss = fnval(f,x(idx));
