function y = sh(n,m,theta,phi);

% y = sh(n,m,theta,phi);
% spherical harmincs Ynm(theta,phi)



a=sqrt( ((2*n+1)/(4*pi)) * factorial(n-m) / factorial(n+m) );

Pvector = legendre(n,cos(theta));
if m>0,
    P = Pvector(m+1,:);
else
    P = Pvector(abs(m)+1,:) * (-1)^abs(m) * factorial(n-abs(m)) / factorial(n+abs(m));
end;

y = a * P * exp(j*m*phi);


    