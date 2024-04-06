function R = random_rotation_matrix()

% Author: Tom Shlomo, ACLab BGU, 2020

x = rand_on_sphere(1, "cart");
theta = rand()*2*pi;
R = axang2rotm([x theta]);

end

