function [th,ph,a] = gaussian(N)
%GAUSSIAN_SAMPLING generates Gaussian distribution of samples 
% on a sphere.
% [a,th,ph]=gaussian_sampling(N);
% N is the order.
% a are the sampling weights.
% th are the elevation angles for all samples.
% ph are the azimuth angles for all samples.
%
% Fundmentals of Spherical Array Processing
% Boaz Rafaely, 2017.

% Find roots of Legendre polynomial
x = roots(legendre_coefficients(N+1));
th = acos(x)';
th = sort(th);

% Sampling weights
a = pi/(N+1) * 2/(N+2)^2 * (1-cos(th).^2) ./ polyval(legendre_coefficients(N+2),cos(th)).^2;

% Equal-angle over azimuth
ph = 0:pi/(N+1):2*pi*(1-1/(2*N+2));

% Build vectors
th=repmat(th,2*N+2,1);
th=reshape(th,1,(N+1)*(2*N+2));
a=repmat(a,2*N+2,1);
a=reshape(a,1,(N+1)*(2*N+2));
ph=repmat(ph,1,N+1);

if nargout==1
    th = [th ph];
end


