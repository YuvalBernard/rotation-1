function [t_ss, Mz_ss] = blochSS(t,Mz)
% This function gets the z-magnetization and calculates the magnetization
% in steady state and the time needed to reach it.
% Fit spline to data and calculate derivative.
f = spline(t,Mz);
dfdt = fnder(f);
% evaluate derivatives at all t
dMzdt = fnval(dfdt,t);
% Find t such that (d/dt)Mz -> 0, arbitrarily (d/dt)Mz < 1e-4.
% Find the 100 first occurences that satisfy condition.
poss_t_ss = find(abs(dMzdt) < 1e-4,100,'first');
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
% Evaluate steady-state z-magnetization based on calculated t_ss.
Mz_ss = fnval(f,t_ss);