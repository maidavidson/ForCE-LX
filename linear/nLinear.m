function [n]=nLinear(k,d)
% [n]=nLinear(k,d)
% Function to compute liner wave theory n-value
% Mark Davidson
kd2=2*k.*d;
n=0.5*(1+(kd2./sinh(kd2)));