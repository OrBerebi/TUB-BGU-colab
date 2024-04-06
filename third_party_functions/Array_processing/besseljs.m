function y = besseljs(n,x);

% spherical Bessel function of the first kind (j)
% computed from 
% Bessel function of the first kind (J)
%
% y - output
% n - order
% x - input

% avoid division by zero
eps=1e-20;
x = x + (x==0).*eps;

% compute output
y = sqrt(pi./(2*x)) .* besselj(n+0.5,x);