function Ygrad = shgrad(omega, N, isComplex, Y)
arguments
    omega       (:,2) double
    N           (1,1) {mustBeNonnegative}
    isComplex   (1,1) logical = true
    Y                 double = shmat(N, omega, isComplex, false)
end
% Author: Tom Shlomo, ACLab BGU, 2020


G = size(omega,1);
Q = (N+1)^2;
assert(size(Y,2)==Q);
assert(G==size(Y,1));

i = (1:Q);
Y = permute(Y,[3 2 1]);
omega = permute(omega, [2 3 1]);
Ygrad = zeros(2, Q, G);
[n,m] = i2nm(i, isComplex);
if isComplex   
    % phi derivative
    Ygrad(2,:,:) = 1i*m.*Y; 
    
    % theta derivative
    [j, I] = nm2i(n, m+1, isComplex);
    j(I) = 1;
    Ynmp1 = Y(1,j,:);
    Ynmp1(1,I,:) = 0;
    Ygrad(1,:,:) = m.*cot(omega(1,:,:)).*Y + ...
                 sqrt( (n-m).*(n+m+1) ) .* exp(-1i*omega(2,:,:)) .* Ynmp1;
%     Ygrad(1,:,abs(mod2(omega(1,:,:),pi))<=eps(pi)) = 0;
else
    error("Real spherical harmonics are not yet supported");
end

end

