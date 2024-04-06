function [a,th,ph]=equiangle_sampling(N);
%EQUIANGLE_SAMPLING generates Equal-angle distribution of samples 
% on a sphere.
% [a,th,ph]=equiangle_sampling(N);
% N is the order.
% a are the sampling weights.
% th are the elevation angles for all samples.
% ph are the azimuth angles for all samples.
%
% Fundmentals of Spherical Array Processing
% Boaz Rafaely, 2017.

L=2*(N+1); % number of samples over each angle

theta = (0:pi/L:(L-1)*pi/L)+(pi/(2*L)); % symmetric along theta
phi   = 0:2*pi/L:(L-1)*2*pi/L;
th=reshape(repmat(theta,L,1),[1,L^2]);
ph=repmat(phi,1,L);

q=(0:N)';
S0=sin((2*q+1)*th);
S=(1./(2*q'+1))*S0;

a = ((8*pi)/(L^2)) * S .*sin(th);

