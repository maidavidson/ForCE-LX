function [Hd]=HoLinear(H1,T,d1,d);
% [Hd]=HoLinear(H,T,d1,d);
% Hd - wave height at depth d, H1 - Reference wave height at depth d1, T - Period, d - depth
% where deep water wave height is required (enter inf for deep water)
% NB Assumes normal incidence over planar bed
% Mark Davidson
if d ==inf; d=10.^17;  end
if d1==inf; d1=10.^17; end
c1 =cLinear(T,d1);
cd =cLinear(T,d);
k1 =kLinear(T,d1);
kd =kLinear(T,d);
n1 =nLinear(k1,d1);
nd =nLinear(kd,d);
Hd=sqrt((n1.*c1) ./ (nd.*cd)) .* H1;
