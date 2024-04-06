function D=wignerD(N,alpha,beta,gamma);
%WIGNERD returns the Wigner-D matrix.
% D = WignerD(N,alpha,beta,gamma);
% N is the order.
% alpha, beta, gamma are the rotation angles.
% Comment: this is a simple but not the most efficient implementation.
%
% Fundmentals of Spherical Array Processing
% Boaz Rafaely, 2017.

D=zeros((N+1)^2,(N+1)^2);

for n=0:N,
    for m=-n:n,
        for mm=-n:n,
            mu=abs(mm-m);
            ni=abs(mm+m);
            s=n-(mu+ni)/2;
            xi=((-1)^(m-mm))^(m<mm);
            c=sqrt(factorial(s)*factorial(s+mu+ni)/factorial(s+mu)/factorial(s+ni));
            d = xi * c * sin(beta/2)^mu * cos(beta/2)^ni * jacobiP(s,mu,ni,cos(beta));
            Dx = exp(-j*mm*alpha) * d * exp(-j*m*gamma);
            D(n^2+n+mm+1,n^2+n+m+1) = Dx;
        end;
    end;
end;

