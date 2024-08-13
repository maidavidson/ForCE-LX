function [c]=cLinear(T,d);
% [c]=cLinear(T,d);
% Approximate solution to the linear wave dispersion equation after Hunt 1979
% Mark Davidson
if isinf(d); d=8000; end
g=9.81; % Gravity
k=kLinear(T,d); % Wave number
c=(g*T./(2*pi)).*tanh(k.*d);