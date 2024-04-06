function y = besselhs(n,x);

% spherical Bessel function of the third kind (h, or spherical Henkel function)
% computed from 
% Bessel function of the third kind (H, or Henkal function)
%
% y - output
% n - order
% x - input

y = sqrt(pi./(2*x)) .* besselh(n+0.5,x);