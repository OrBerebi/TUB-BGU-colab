function t = angle_between(omega1, omega2)
arguments
    omega1 (:,2) double
    omega2 (:,2) double
end
% Author: Tom Shlomo, ACLab BGU, 2020


t = real(acos( s2c(omega1)*s2c(omega2)' ));
% Equivalent mathematically to, but more numerically stable than:
% t = acos( cos(omega1(:,1).*cos(omega2(:,1)) + cos(omega1(:,2)-omega2(:,2)).*sin(omega1(:,1)).*sin(omega2(:,1))) );

end

