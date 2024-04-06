function y = besselhsd(n,x);
%BESSELHSD return the first derivative of the spherical Hankel function 
% of the first kind.
% y = besselhsd(n,x);
% n is the order.
% x is the argument.
%
% Fundmentals of Spherical Array Processing
% Boaz Rafaely, 2017.

% Avoid division by zero
x = x + (x==0).* 1e-20;

y = (n./x) .* besselhs(n,x) - besselhs(n+1,x);