function y = besseljsd(n,x);

% First derivative of the
% spherical Bessel function of the first kind (j)
% computed from 
% spherical Bessel function of the first kind
%
% y - output
% n - order
% x - input

y = (n./x) .* besseljs(n,x) - besseljs(n+1,x);