function [E]=ELinear(H);
% [E]=ELinear(H);
% Linear Wave Energy Function
% H - wave height, E- Wave energy 
% Mark Davidson
E=(1025*9.81*H.^2)./16;