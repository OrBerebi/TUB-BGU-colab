function Y = sh2(N,theta,phi);
%SH2 returns the spherical harmonics matrix.
% Y = sh2(N,theta,phi);
% N is the maximum order.
% theta are the elevation angles of all points.
% phi are the azimuth angles of all points.
% Y is the spherical harmonics matrix, size (N+1)^2 by length(theta).
%
% Fundmentals of Spherical Array Processing
% Boaz Rafaely, 2017.

% complex constant
j=sqrt(-1);

% Number of points on the sphere
L=length(theta);

% Generate row vectors
theta=reshape(theta,1,L); 
phi=reshape(phi,1,L);     

% n=0
Y=sqrt(1/(4*pi))*ones(1,L);

for n=1:N,
    
    % positive m
    for m=0:n,
        a(m+1) = sqrt( ((2*n+1)/(4*pi)) * factorial(n-m) / factorial(n+m) );
    end;
    a=reshape(a,n+1,1);
    m=[0:n]';
    Y1 = (a*ones(1,L)) .* legendre(n,cos(theta)) .* exp(j*m*phi) ;
        
    % negative m
    m=[-n:-1]';
    Y2 = (((-1).^m)*ones(1,L)) .* conj(Y1(end:-1:2,:));
            
    % append to Y
    Y = [Y; Y2; Y1];
    
end;
