function T = legendre_coefficients(N);
%LEGENDRE_COEFFICIENTS returns the Legendre polynomial coefficients.
% T = legendre_coefficients(N);
% N is the order. 
% The coefficients are computed using the binomial formula, see, for example, 
% Arfken et al, Mathematical Methods for Physicists.
%
% Fundmentals of Spherical Array Processing
% Boaz Rafaely, 2018.

T=zeros(1,N+1);
for r=0:floor(N/2),
    T(2*r+1)=(1/2^N) * (-1)^r * factorial(2*N-2*r) / (factorial(r)*factorial(N-r)*factorial(N-2*r) );
end;


