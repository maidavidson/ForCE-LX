function [Lo]=LoLinear(T);
% deep water wavelength
% Mark Davidson
g=9.81;     % gravity
Lo=g*T.^2/(2*pi);




