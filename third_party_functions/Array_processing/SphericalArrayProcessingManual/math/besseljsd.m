function y = besseljsd(n,x);
%BESSELJSD returns the first derivative of the spherical Bessel function 
% of the first kind.
% y = besseljsd(n,x);
% n is the order.
% x is the argument.
%
% Fundmentals of Spherical Array Processing
% Boaz Rafaely, 2017.

% Avoid division by zero
x = x + (x==0).* 1e-20;

y = (n./x) .* besseljs(n,x) - besseljs(n+1,x);