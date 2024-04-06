function Y = sh2(N,theta,phi)
% function Y = sh2(N,theta,phi);
%
% N maximum order for n
% theta are the angles theta for the entire location vector
% phi are the angle phi for the entire location vector
% Y is (N+1)^2 by length(theta)

% complex constant
j=sqrt(-1);

% ensure sizes are correct
if length(theta) ~= length(phi);
    fprintf('Lengths of theta and phi must be equal!');
    return;
end;
L=numel(theta);
theta=reshape(theta,1,L); % make a row vector
phi=reshape(phi,1,L);       % make a row vector


%n=0
Y=sqrt(1/(4*pi))*ones(1,L);

for n=1:N,
    
    % positive m
    a=zeros(1,n);
    for m=0:n,
        a(m+1) = sqrt( ((2*n+1)/(4*pi)) * factorial(n-m) / factorial(n+m) );
    end;
    a=reshape(a,n+1,1);
    m=(0:n)';
    Y1 = (a*ones(1,L)) .* legendre(n,cos(theta)) .* exp(j*m*phi) ;
        
    % negative m
    m=(-n:-1)';
    Y2 = (((-1).^m)*ones(1,L)) .* conj(Y1(end:-1:2,:));
            
    % append to Y
    Y = [Y; Y2; Y1];
    
end;
