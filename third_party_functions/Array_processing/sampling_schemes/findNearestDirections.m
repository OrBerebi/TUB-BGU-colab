%findNearestDirections - finds the nearest directions from a given directions
%% INPUTS:
%  hrtf_directions            : data for direction samples, array of dimension = num_directions x (2 for azim, elev)
%  given_directions           : data for given direction, array of dimension = num_directions x (2 for azim, elev)
%  visualize                  : flag for figure
%
%% OUTPUTS:
%  NEAREST_DIRS               : the nearest found directions from the original hrtf_directions
%  index                      : return the indexs of the original taken directions
function [ NEAREST_DIRS, indexes, errors] = findNearestDirections(hrtf_directions, given_directions, visualize)

num_directions = size(given_directions,1);
num_HRTF_dirs = size(hrtf_directions,1); 

hrtf_directions(:,2) = pi/2 - hrtf_directions(:,2);
given_directions(:,2) = pi/2 - given_directions(:,2);

if size(hrtf_directions,2)<3
    hrtf_directions(:,3) = ones(size(hrtf_directions,1),1);
end
if size(given_directions,2)<3
    given_directions(:,3) = ones(size(given_directions,1),1);
end

[hrtf_dir_x, hrtf_dir_y, hrtf_dir_z] = sph2cart(hrtf_directions(:,1),hrtf_directions(:,2),hrtf_directions(:,3));
hrtf_dirs_cart = [hrtf_dir_x hrtf_dir_y hrtf_dir_z];

NEAREST_DIRS = [];

for i=1:1:num_directions
    
    given_dir = given_directions(i,:);
    
    % dot product metric
    [given_dir_x, given_dir_y, given_dir_z] = sph2cart(given_dir(1), given_dir(2), given_dir(3));
    given_dir_cart = [given_dir_x, given_dir_y, given_dir_z];
    given_dirs_cart = ones(num_HRTF_dirs,1) * given_dir_cart;
    dott = dot(hrtf_dirs_cart, given_dirs_cart, 2);
    [~, indexes(i)] = max(dott);
    errors(i) = norm(given_dir_cart - hrtf_dirs_cart(indexes(i),:));
    % nearest direction in the sampled HRTF dataset
    nearest_dir = hrtf_directions(indexes(i),:);
    NEAREST_DIRS = [NEAREST_DIRS; nearest_dir];
    
end

% visualize HRTF samples, given direction, nearest direction
if visualize
    SIZE = 10 * ones(num_directions,1);
    SIZE_hrtf = 2 * ones(num_HRTF_dirs,1);
    [X,Y,Z] = sph2cart(NEAREST_DIRS(:,1),NEAREST_DIRS(:,2),NEAREST_DIRS(:,3));
    [X2,Y2,Z2] = sph2cart(given_directions(:,1),given_directions(:,2),given_directions(:,3));
    figure;
    scatter3(hrtf_dir_x, hrtf_dir_y, hrtf_dir_z, SIZE_hrtf);
    hold on
    scatter3(X2, Y2, Z2, SIZE*2, 'filled');
    scatter3(X, Y, Z, SIZE*2, 'filled');
    quiver3(0,0,0,1,0,0, 1.75, 'filled', 'r', 'LineWidth', 2);   % X
    quiver3(0,0,0,0,1,0, 1.75, 'filled', 'g', 'LineWidth', 2);   % Y
    quiver3(0,0,0,0,0,1, 1.75, 'filled', 'b', 'LineWidth', 2);   % Z
    xlabel('X');    ylabel('Y');    zlabel('Z');
    axis equal;
    view([40 20]),    %title(label);
    legend('Original HRTF directions','Given directions','Found directions','Location','north')
end

NEAREST_DIRS(:,2) = pi/2- NEAREST_DIRS(:,2);
end

