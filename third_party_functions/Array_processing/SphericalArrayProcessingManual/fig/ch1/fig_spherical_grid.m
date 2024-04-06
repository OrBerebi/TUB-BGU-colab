%FIG_SPHERICAL_GRID generates Figure 1.9 presenting a sphere with an 
% angular grid.
%
% Fundmentals of Spherical Array Processing
% Boaz Rafaely, 2017.

close all;
clear all;

figure(1);
[x,y,z] = sphere(16);
h1=surf(x,y,z);
set(h1,'EdgeColor',[0 0 0.5]); % dark blue
set(h1,'LineWidth',1.5);
colormap(white);
axis tight equal 
axis off
set(h1,'EdgeAlpha',1);
set(h1,'FaceAlpha',0.8);
zoom(1);

% print -dpng -r400 ../../figures/chapter01/fig_sphere_grid.png
