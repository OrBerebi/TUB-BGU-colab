function out = shMatrix(N,el,az)
% function out = shMatrix(N,el,az)
%
% Generate a matrix of spherical harmonics
%
%    N   - Maximum order
%    el  - Elevation vector
%    az  - Azimuth vector
% 
%    out - (N+1)^2 x length(el)  matrix
%
%  July 2014, Boaz Rafaely, Ben-Gurion University
%  edited by Jonathan Sheaffer
%  Part of the EARS Beamforming Toolbox
%
% complex constant
j=sqrt(-1);

% ensure sizes are correct
if length(el) ~= length(az);
    fprintf('Lengths of az and el must be equal!');
    return;
end;
L=numel(el);
theta=reshape(el,1,L); % make a row vector
phi=reshape(az,1,L);       % make a row vector


%n=0
out=sqrt(1/(4*pi))*ones(1,L);

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
    out = [out; Y2; Y1];
    
end