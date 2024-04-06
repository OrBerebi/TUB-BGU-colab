function [omega, disp] = mean_doa(omega, weights)
arguments
    omega (:,2) double
    weights (:,1) double {mustBeNonnegative} = ones(size(omega,1),1)
end
% Author: Tom Shlomo, ACLab BGU, 2020


x = s2c(omega);
x = sum(x.*weights,1);
omega = c2s(x);
omega = omega(1:2);

if nargout>=2
    disp = 1 - x(3)/sum(weights);
end

end

