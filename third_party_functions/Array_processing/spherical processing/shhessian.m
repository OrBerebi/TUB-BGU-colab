function Yhessian = shhessian(omega, N, isComplex, Y, Ygrad)
arguments
    omega       (:,2)   double
    N           (1,1)   {mustBeNonnegative}
    isComplex   (1,1)   logical = true
    Y                   double = shmat(N, omega, isComplex, false)
    Ygrad       (2,:,:) double = shgrad(omega, N, isComplex, Y)
end
% Author: Tom Shlomo, ACLab BGU, 2020


G = size(omega,1);
Q = (N+1)^2;
assert(size(Y,2)==Q);
assert(G==size(Y,1));
assert(Q==size(Ygrad,2));
assert(G==size(Ygrad,3));

i = permute((1:Q), [1 3 2]);        % 1 x 1 x Q x 1
Y = permute(Y,[3 4 2 1]);           % 1 x 1 x Q x G
Ygrad = permute(Ygrad, [1 4 2 3]);  % 2 x 1 x Q x G
omega = permute(omega, [2 3 4 1]);  % 2 x 1 x 1 x G
Yhessian = zeros(2, 2, Q, G);       % 2 x 2 x Q x G
[n,m] = i2nm(i, isComplex);         % 1 x 1 x Q x 1
if isComplex   
    % phi phi derivative
    Yhessian(2,2,:,:) = -m.^2 .* Y; 
    % phi theta derivative
    [j, I] = nm2i(n, m+1, isComplex);
    j(I) = 1;
    Ynmp1 = Y(1,1,j,:);
    Ynmp1(1,1,I,:) = 0;
    Ygradnmp1 = Ygrad(:, :, j, :);
    Ygradnmp1(:, :, I, :) = 0;
    C = sqrt( (n-m).*(n+m+1) ) .* exp(-1i*omega(2,:,:,:));
    Yhessian(2,1,:,:) = 1i* m   .*( m.*cot(omega(1,:,:,:)).*Y ) + ...
                        C.*( Ygradnmp1(2,1,:,:) + -1i*Ynmp1 );
    Yhessian(1,2,:,:) = Yhessian(2,1,:,:);
    
    % theta theta derivative
    Yhessian(1,1,:,:) = m.*(  csc(omega(1,:,:,:)).^2  .* Y ...
                            + cot(omega(1,:,:,:))     .* Ygrad(1, 1, :, :) ) + ...
                        C.*Ygradnmp1(1, 1, :, :);
else
    error("Real spherical harmonics are not yet supported");
end

end

