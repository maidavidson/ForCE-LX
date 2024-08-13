function [k]=kLinear(T,d);
% [k]=kLinear(T,d);
% Approximate solution to the linear wave dispersion equation after Hunt 1979
% Mark Davidson
g=9.81;     % gravity
y=4.03*d./(T.^2);
kd2=y.^2 + y ./ ...
    (1+(0.666.*y)+(0.355.*y.^2)+(0.161.*y.^3)+(0.0632.*y.^4)+(0.0218.*y.^5)+(0.00564.*y.^6));
k=sqrt(kd2)./d;



