function [g, omega] = sphere_abs_newton(fnm, omega0, NameValueArgs)
arguments
    fnm (:,1) double
    omega0 (1,2) double
    NameValueArgs.tol (1,1) double {mustBePositive, mustBeFinite} = 2*pi*1e-6
    NameValueArgs.maxIter (1,1) double {mustBeInteger, mustBePositive} = 50
    NameValueArgs.isComplex (1,1) logical = true
    NameValueArgs.plotFlag (1,1) logical = false
    NameValueArgs.parent (1,1) matlab.graphics.axis.Axes
    NameValueArgs.Y0 double = [];
    NameValueArgs.f0 double = [];
    NameValueArgs.stopIfNotNegativeDefinite (1,1) logical = false
    NameValueArgs.stopIfNotPositiveDefinite (1,1) logical = false
end
% Author: Tom Shlomo, ACLab BGU, 2020


N = sqrt(size(fnm,1))-1;
omega = nan(NameValueArgs.maxIter+1, 2);
omega(1,:) = omega0;
g = nan(NameValueArgs.maxIter,1);

for i=1:NameValueArgs.maxIter
    if i==1 && ~isempty(NameValueArgs.Y0)
        Y = NameValueArgs.Y0;
    else
        Y = shmat(N, omega(i,:), NameValueArgs.isComplex, false);
    end
    if i==1 && ~isempty(NameValueArgs.f0)
        f = NameValueArgs.f0;
    else
        f = Y*fnm;
    end
    grad_Y = shgrad(omega(i,:), N, NameValueArgs.isComplex, Y);
    grad_f = grad_Y*fnm;
    H_Y = shhessian(omega(i,:), N, NameValueArgs.isComplex, Y, grad_Y);
    H_f = sum( H_Y .* reshape(fnm, 1, 1, []), 3);
    
    g(i) = abssq(f);
    grad_g = 2*( real(f)*real(grad_f) + imag(f)*imag(grad_f) );
    H_g = 2*( real(grad_f)*real(grad_f)' + real(f)*real(H_f) + ...
              imag(grad_f)*imag(grad_f)' + imag(f)*imag(H_f) );
   
    grad_z = grad_g / g(i);
    H_z = H_g / g(i) - ( grad_g * grad_g' / g(i)^2 );
    if NameValueArgs.stopIfNotNegativeDefinite  && ~(trace(H_z)<0 && det(H_z)>0)
        break
    end
    if NameValueArgs.stopIfNotPositiveDefinite  && ~(trace(H_z)>0 && det(H_z)>0)
        break
    end
    domega = - H_z \ grad_z;
    omega(i+1,:) = omega(i,:) + domega.';
    if angle_between(omega(i,:), omega(i+1,:)) <= NameValueArgs.tol
        i = i+1; %#ok<FXSET>
        Y = shmat(N, omega(i,:), NameValueArgs.isComplex, false);
        f = Y*fnm;
        g(i) = abssq(f);
        break
    end
end
omega = omega(1:i,:);
g = g(1:i,:);

if NameValueArgs.plotFlag
    if ~isfield(NameValueArgs, 'parent')
        NameValueArgs.parent = gca;
    end
    hammer.surf([],fnm, [], @(f) 10*log10(abssq(f)), true, 'Parent', NameValueArgs.parent);
    hold on;
    hammer.plot(omega(:,1), omega(:,2), 'r.-', 'Parent', NameValueArgs.parent);
    hold off;
end
end

