function [P]=PLinear(E,C,n);
% [P]=PLinear(E,C,n);
% Linear Wave Power Function
% E - Wave energy, C - wave celerity, n - linear n-value 
% Mark Davidson
P=E.*C.*n;