function y = besselhs(n,x);
%BESSELHS return the spherical Hankel function of the first kind.
% y = besselhs(n,x);
% n is the order.
% x is the argument.
%
% Fundmentals of Spherical Array Processing
% Boaz Rafaely, 2017.

% Avoid division by zero
x = x + (x==0).* 1e-20;

y = sqrt(pi./(2*x)) .* besselh(n+0.5,x);