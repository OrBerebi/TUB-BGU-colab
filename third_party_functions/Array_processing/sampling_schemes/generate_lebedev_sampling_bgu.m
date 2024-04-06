function [a,th,ph,max_SHorder] = generate_lebedev_sampling_bgu(num_points, minTheta, maxTheta, plot)

% [a,th,ph,max_SHorder] = generate_lebedev_sampling(num_points, minTheta, maxTheta, plot)
% 
% Lebedev sampling scheme
% 
% Inputs:
%         num_points :          wanted number of points (the grid will be as close as possible to this number
%         minTheta   :          minimum elevation for the HRTF measurement (in degree, between -90 to 90)
%         maxTheta   :          maximum elevation for the HRTF measurement (in degree, between -90 to 90)
%         plot       :          flag for plotting the sampling point
%
% Outputs:
%         a          :          weights for discrete integration
%         th         :          vector of elevation angles (in radians, -pi/2 to pi/2)
%         ph         :          vector of azimuth angles   (in radians, 0 to 2pi)
%         max_SHorder:          maximun spherical harmonics order for the scheme
%
% External routines needed: sofia_lebedev.m, lebedev_calc.m
%

if nargin<4
    plot = 0; % default dont plot
end
degrees_avail=[6, 14, 26, 38, 50, 74, 86, 110, 146, 170, 194, 230,...
    266, 302, 350, 434, 590, 770, 974, 1202, 1454, 1730, 2030, 2354, 2702,...
    3074, 3470, 3890, 4334, 4802, 5294, 5810];

[~,~,V] = findnearest(num_points,degrees_avail);
[gridData, ~, max_SHorder] = sofia_lebedev(V(1), 0);

a = (gridData(:,3))*4*pi;
th = pi/2 - gridData(:,2);
ph = gridData(:,1);

ind_th = find(th>=minTheta*pi/180 & th<=maxTheta*pi/180);

th = th(ind_th);
ph = ph(ind_th);
a = a(ind_th);
if plot    
    plotSampleDirections([ph th],[],' ',1,1,8); colorbar off
    colormap Gray;
    hold on;
    sphere;
    axis equal;
    view([40 -10]); 
    rotate3d on;
    light;
    alpha(.8);
    lighting phong;
    hold off;
end