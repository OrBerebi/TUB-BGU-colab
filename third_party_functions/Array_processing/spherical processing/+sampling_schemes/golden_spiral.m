function [th, ph] = golden_spiral(K, applyRandomRotation)
arguments
    K (1,1) double {mustBePositive, mustBeFinite, mustBeInteger}
    applyRandomRotation (1,1) logical = false;
end
% Author: Tom Shlomo, ACLab BGU, 2020

%%
i = (0:(K-1))' + 0.5;
th = acos(1-2*i/K);
ph = pi*(1+sqrt(5))*i;

%%
if applyRandomRotation
    x = s2c(th, ph);
    R = random_rotation_matrix();
    x = x*R;
    [th, ph] = c2s(x(:,1), x(:,2), x(:,3));
end

%%
if nargout==1
    th = [th ph];
end

end