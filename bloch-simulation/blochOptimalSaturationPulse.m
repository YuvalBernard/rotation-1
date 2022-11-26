function [opt_w1, opt_dw, opt_t_ss, opt_Mz_ss, T] = blochOptimalSaturationPulse(t,T1,T2,w1,dw)
% This function finds the optimal saturation pulse
% Saturation pulse is defined such that at t_ss, Mz_ss -> 0.
% We want to find the shortest saturation pulse obtainable.
% We can define Mz_ss->0 such that Mz<=0.015. Find the indices of Mz_ss
% that satisfy the condition.
% Then check for which of these indices t_ss is minimal.
% Get values of w1 and dw that satisfy the condition.

% Calculate t_ss and Mz_ss for all given combinations of dw,w1.
[DW, W1] = meshgrid(dw,w1);
[t_ss,Mz_ss] = arrayfun(@(w1,dw) blochSS(t,bloch(t,T1,T2,w1,dw)), W1, DW);

% Find which Mz_ss are saturated.
idx_sat = find(Mz_ss <= 0.015);

% Check for which of these values t_ss is minimal.
sel_t_ss = t_ss(idx_sat);
[opt_t_ss, idx_opt] = min(sel_t_ss);
% Before we return opt_w1 and opt_dw, check first if sel_t_ss contains
% duplicate values of opt_t_ss. In case of duplictes, idx_opt shall be the
% one for which w1 is minimal (conserve more power).

% Get positions of duplicates (if present)
poss_opt_t_ss = sel_t_ss == opt_t_ss;
% If no duplicates, conclude:
if length(idx_sat(poss_opt_t_ss)) == 1
    opt_w1 = W1(idx_sat(idx_opt));
    opt_dw = DW(idx_sat(idx_opt));
    opt_Mz_ss = Mz_ss(idx_sat(idx_opt));
else
% Find minimal w1 in subset of w1 that satisfies opt_t_ss
    sel_w1 = W1(idx_sat(poss_opt_t_ss));
    [opt_w1, idx_opt] = min(sel_w1);
    opt_dw = DW(idx_sat(idx_opt));
    opt_Mz_ss = Mz_ss(idx_sat(idx_opt));
end

% Also return a table containing all values of dw,w1
% for which Mz_ss is saturated and also dw != 0. (because that's trivial.)
alt_idx = idx_sat(DW(idx_sat ~= 0));
T = table([opt_w1;W1(alt_idx)],[opt_dw;DW(alt_idx)],...
    [opt_t_ss;t_ss(alt_idx)],[opt_Mz_ss;Mz_ss(alt_idx)],...
    'VariableNames',{'w_1','dw','t^{ss}','Mz^{ss}'});
