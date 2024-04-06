function vnm_th = derivative_th(vnm,th,ph);
%DERIVATIVE_TH returns the derivative along theta of a spherical harmonics
% steering vector.
% vnm_th = derivative_th(vnm,th,ph);
% vnm is the steering vector.
% theta and phi represent position point of the derivative on the sphere.
%
% Fundmentals of Spherical Array Processing
% Boaz Rafaely, 2017.

N=sqrt(length(vnm))-1;

for n=0:N, 
    for m=-n:n,
        g1(n^2+n+m+1,1)=m*cot(th);
        g2(n^2+n+m+1,1)=sqrt((n-m)*(n+m+1))*exp(j*ph);
    end;
end

for n=0:N,
    for m=-n:n,
        q=n^2+n+m+1;
        vnm_th(q,1)=g1(q)*vnm(q);
        if m<n, vnm_th(q,1)=vnm_th(q,1)+g2(q)*vnm(q+1); end
    end;
end;

