function [ NEAREST_DIRS , given_directions, index] = findPseudoUniformSamplingPoints_bgu(hrtf_directions, N ,visualize)
% findPseudoUniformSamplingPoints
% The functions find the nearest points from a given distribution over the
% sphere to a uniform distribution with order N
%
% INPUTS:
%  hrtf_directions            : data for direction samples, array of dimension = num_directions x (2 for azim [0 2pi], elev[0 pi])
%  N                          : wanted SH order for the uniform sampling (avaliable only for N = 2-10)
%  visualize                  : flag for figure
%
% OUTPUTS:
%  NEAREST_DIRS               : the nearest found directions from the original hrtf_directions
%  given_directions           : the uniform sampling directions
%  index                      : return the indexs of the original taken directions

[~,el,az]=uniform_sampling_extended(N);
given_directions = [az.'  el.'];

num_directions = length(el);
num_HRTF_dirs = size(hrtf_directions,1);

[hrtf_dir_x, hrtf_dir_y, hrtf_dir_z] = sph2cart(hrtf_directions(:,1),pi/2 - hrtf_directions(:,2),1);
hrtf_dirs_cart = [hrtf_dir_x hrtf_dir_y hrtf_dir_z];

NEAREST_DIRS = [];

for i=1:1:num_directions
    
    given_dir = given_directions(i,:);
    
    % dot product metric
    [given_dir_x, given_dir_y, given_dir_z] = sph2cart(given_dir(1), pi/2 - given_dir(2),1);
    given_dir_cart = [given_dir_x, given_dir_y, given_dir_z];
    given_dirs_cart = ones(num_HRTF_dirs,1) * given_dir_cart;
    dott = dot(hrtf_dirs_cart, given_dirs_cart, 2);
    [~, index(i)] = max(dott);
    
    % nearest direction in the sampled HRTF dataset
    nearest_dir = hrtf_directions(index(i),:);
    NEAREST_DIRS = [NEAREST_DIRS; nearest_dir];
    
end

% visualize HRTF samples, given direction, nearest direction
if visualize
    SIZE = 10 * ones(num_directions,1);
    SIZE_hrtf = 2 * ones(num_HRTF_dirs,1);
    [X,Y,Z] = sph2cart(NEAREST_DIRS(:,1),pi/2 - NEAREST_DIRS(:,2),1);
    [X2,Y2,Z2] = sph2cart(given_directions(:,1),pi/2 - given_directions(:,2),1);
    figure;
    scatter3(hrtf_dir_x, hrtf_dir_y, hrtf_dir_z, SIZE_hrtf);
    hold on
    scatter3(X2, Y2, Z2, SIZE*2, 'filled');
    scatter3(X, Y, Z, SIZE*2, 'filled');
    quiver3(0,0,0,1,0,0, 1.75, 'filled', 'r', 'LineWidth', 2);   % X
    quiver3(0,0,0,0,1,0, 1.75, 'filled', 'g', 'LineWidth', 2);   % Y
    quiver3(0,0,0,0,0,1, 1.75, 'filled', 'b', 'LineWidth', 2);   % Z
    xlabel('X');    ylabel('Y');    zlabel('Z');
    title(['Nearest uniform points of order N=',num2str(N), '. Total points ',num2str(num_directions)])
    axis equal;
    view([40 20]),    %title(label);
    legend('Original HRTF directions','Uniform sampling','Found directions','Location','north')
end


end

