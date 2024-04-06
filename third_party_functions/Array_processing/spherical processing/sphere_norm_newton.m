function [g, omega, notDefiniteFlag] = sphere_norm_newton(Fnm, omega0, NameValueArgs)
arguments
    Fnm (:,:) double
    omega0 (1,2) double
    NameValueArgs.tol (1,1) double {mustBePositive, mustBeFinite} = 2*pi*1e-6
    NameValueArgs.maxIter (1,1) double {mustBeInteger, mustBePositive} = 50
    NameValueArgs.isComplex (1,1) logical = true
    NameValueArgs.plotFlag (1,1) logical = false
    NameValueArgs.parent (1,1) matlab.graphics.axis.Axes = gca
    NameValueArgs.Y0 double = [];
    NameValueArgs.F0 double = [];
    NameValueArgs.stopIfNotNegativeDefinite (1,1) logical = false
    NameValueArgs.stopIfNotPositiveDefinite (1,1) logical = false
end
% Author: Tom Shlomo, ACLab BGU, 2020


N = ceil(sqrt(size(Fnm,1))-1);
K = size(Fnm,2);
omega = nan(NameValueArgs.maxIter+1, 2);
omega(1,:) = omega0;
g = nan(NameValueArgs.maxIter,1);
notDefiniteFlag = false;
for i=1:NameValueArgs.maxIter
    if i==1 && ~isempty(NameValueArgs.Y0)
        Y = NameValueArgs.Y0;
    else
        Y = shmat(N, omega(i,:), NameValueArgs.isComplex, false);
    end
    if i==1 && ~isempty(NameValueArgs.F0)
        F = NameValueArgs.F0;
    else
        F = Y*Fnm; % Q x 1
    end
    grad_Y = shgrad(omega(i,:), N, NameValueArgs.isComplex, Y);
    grad_F = grad_Y*Fnm; % 2 x Q
    H_Y = shhessian(omega(i,:), N, NameValueArgs.isComplex, Y, grad_Y); % 2 x 2 x sh
    g(i) = F*F';
    grad_g = zeros(2,1);
    H_g = zeros(2);
    for k=1:K
        fnm = Fnm(:,k);
        f = F(k);
        grad_f = grad_F(:,k);
        H_f = sum( H_Y .* reshape(fnm, 1, 1, []), 3);
        grad_g = grad_g + 2*( real(f)*real(grad_f) + imag(f)*imag(grad_f) );
        H_g = H_g + 2*( real(grad_f)*real(grad_f)' + real(f)*real(H_f) + ...
                  imag(grad_f)*imag(grad_f)' + imag(f)*imag(H_f) );
    end
    grad_z = grad_g / g(i);
    H_z = H_g / g(i) - ( grad_g * grad_g' / g(i)^2 );
    if NameValueArgs.stopIfNotNegativeDefinite  && ~(trace(H_z)<0 && det(H_z)>0)
        notDefiniteFlag = true;
        break
    end
    if NameValueArgs.stopIfNotPositiveDefinite  && ~(trace(H_z)>0 && det(H_z)>0)
        notDefiniteFlag = true;
        break
    end
    domega = - H_z \ grad_z;
    omega(i+1,:) = omega(i,:) + domega.';
    if angle_between(omega(i,:), omega(i+1,:)) <= NameValueArgs.tol
        i = i+1; %#ok<FXSET>
        Y = shmat(N, omega(i,:), NameValueArgs.isComplex, false);
        F = Y*Fnm;
        g(i) = F*F';
        break
    end
end
omega = omega(1:i,:);
g = g(1:i,:);

if NameValueArgs.plotFlag
    hammer.surf([],Fnm, [], @(f) f, true, 'Parent', NameValueArgs.parent);
    hold on;
    hammer.plot(omega(:,1), omega(:,2), 'r.-', 'Parent', NameValueArgs.parent);
    hold off;
end
end

