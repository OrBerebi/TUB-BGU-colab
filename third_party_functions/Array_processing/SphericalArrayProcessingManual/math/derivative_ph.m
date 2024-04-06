function vnm_ph = derivative_ph(vnm);
%DERIVATIVE_PH returns the derivative along phi of a spherical harmonics
% steering vector.
% vnm_ph = derivative_ph(vnm);
% vnm is the steering vector.
%
% Fundmentals of Spherical Array Processing
% Boaz Rafaely, 2017.

N=sqrt(length(vnm))-1;
for n=0:N, 
    for m=-n:n, 
        g(n^2+n+m+1,1)=-j*m; 
    end; 
end;

vnm_ph=g.*vnm;
