function [x,y,z] = s2c(theta,phi,r);
%S2C converts spherical to Cartesian coordinates (non Matlab notation!).
% [x,y,z] = s2c(theta,phi,r);
% (x,y,z) is the conventional cartezian coordinates.
% theta is the angle going down from the z-axiz.
% phi is the azimuth angle from the poitive x-axis, towards the positive y-axis.
%
% Fundmentals of Spherical Array Processing
% Boaz Rafaely, 2017.

x = r.*sin(theta).*cos(phi);
y = r.*sin(theta).*sin(phi);
z = r.*cos(theta);

