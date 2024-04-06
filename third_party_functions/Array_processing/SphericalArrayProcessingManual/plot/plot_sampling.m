function plot_sampling(th,ph);
%PLOT_SAMPLING generates contour and balloon plots of the 
% sampling scheme.
% plot_sampling(th,ph);
% (th,ph) are the sampling points in spherical coordinates.
%
% Fundmentals of Spherical Array Processing
% Boaz Rafaely, 2017.

[x,y,z] = s2c(th,ph,1);

AxisFontSize=14;

figure;
n=48;
[X,Y,Z] = sphere(n);
h=surf(X,Y,Z);
colormap(bone);
set(h,'EdgeAlpha',0.1);
axis tight equal
axis off
hold on;
v=plot3(x,y,z,'.');
set(v,'MarkerSize',12,'MarkerEdgeColor','k','MarkerFaceColor','k');

figure;
v=plot((180/pi)*ph,(180/pi)*th,'.','Color',[0 0 0.5]);
set(gca,'FontSize',AxisFontSize);
set(v,'MarkerSize',AxisFontSize);
xlabel('$\phi$ (degrees)','FontSize',AxisFontSize,'Interp','Latex');
ylabel('$\theta$ (degrees)','FontSize',AxisFontSize,'Interp','Latex');
axis([0,360,0,180]);
