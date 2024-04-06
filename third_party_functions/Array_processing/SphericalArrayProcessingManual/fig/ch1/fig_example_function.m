%FIG_EXAMPLE_FUNCTION generates Figures 1.2-1.4, which plot a function on
% the sphere using three ploting methods, namely plot on a sphere, contour
% plot, and balloon plot.
%
% Fundmentals of Spherical Array Processing
% Boaz Rafaely, 2017.

close all;
clear all;

path(path,'../../math');
path(path,'../../plot');

AxisFontSize=14;

N=2;
fnm=zeros((N+1)^2,1);
fnm(5)=2*sqrt(2*pi/15); % (n,m)=(2,-2)
fnm(9)=2*sqrt(2*pi/15); % (n,m)=(2,2)

figure;
plot_sphere(fnm);

figure;
plot_contour(fnm,[],0);

figure;
plot_balloon(fnm,[],0.2);
axis on;
xlabel('$x$','FontSize',AxisFontSize,'Interp','Latex');
ylabel('$y$','FontSize',AxisFontSize,'Interp','Latex');
zlabel('$z$','FontSize',AxisFontSize,'Interp','Latex');
set(gca,'FontSize',AxisFontSize);

% figure(1); print -dpng -r600 ../../figures/chapter01/fig_function_sphere.png
% figure(2); print -dpng -r600 ../../figures/chapter01/fig_function_contour.png
% figure(3); print -dpng -r600 ../../figures/chapter01/fig_function_balloon.png





