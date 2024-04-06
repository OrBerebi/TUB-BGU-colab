function [a,th,ph]=generate_gaussian_sampling(N)

% [a,th,ph]=generate_gaussian_sampling(N);
%
% a - weights
% th - theta angles for all points
% ph - phi angles for all points
% Total no. of points is length(a)=length(th)=length(ph)
%
% Gaussina sampling for order N
%
% Boaz Rafaely 12 October 2006


Lg=N+1;
Leq=2*(N+1);

% Equiangle samples over phi
phi=0:2*pi/Leq:(Leq-1)*2*pi/Leq;


% Gaussian samples and nodes over theta

[x,w]=lgwt(N+1,-1,1);
theta=acos([x; -x(end-1:-1:1)]);
aa=(pi/(N+1))*[w; w(end-1:-1:1)];

count=0;
for j=1:Lg,
    for k=1:Leq,
        count=count+1;
        th(count) = theta(j);
        ph(count) = phi(k);
        a(count)  = aa(j);
    end;
end;



