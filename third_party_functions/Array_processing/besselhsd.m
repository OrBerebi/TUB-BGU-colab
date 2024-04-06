function y = besselhsd(n,x);

% First derivative of the
% spherical Bessel function of the third kind (h, or spherical Henkal function)
% computed from 
% spherical Bessel function of the third kind
%
% y - output
% n - order
% x - input

y = (n./x) .* besselhs(n,x) - besselhs(n+1,x);