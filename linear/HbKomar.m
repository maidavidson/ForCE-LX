% [Hb]=HbKomar(Ho,T);
% Hb = breaking wave height (m) 
% Ho = offshore significant wave height (m)
% T  = peak wave period (s)
% P.D. Komar (1974) equation
function [Hb]=HbKomar(Ho,T);
Hb = 0.39*9.81.^0.2*(T.*Ho.^2).^0.4;