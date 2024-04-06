%FIG_XYZ_COORDINATES generates cartesian axes plots for Figures 1.5-1.8 and 
% figure 1.15.
%
% Fundmentals of Spherical Array Processing
% Boaz Rafaely, 2017.

close all
clear all

FSIZE=70;
FNAME='Helvetica';

figure(1);
line([0 1],[0 0],[0 0],'LineWidth',2,'Color','k');
line([0 0],[0 1],[0 0],'LineWidth',2,'Color','k');
line([0 0],[0 0],[0 1],'LineWidth',2,'Color','k');
axis([0 1 0 1 0 1]);
axis off
axis square
text(1.15,0.07,-0.03,'$x$','FontSize',FSIZE,'FontName',FNAME,'Interp','Latex');
text(-0.02,1.38,-0.13,'$y$','FontSize',FSIZE,'FontName',FNAME,'Interp','Latex');
text(-0.12,0,1.3,'$z$','FontSize',FSIZE,'FontName',FNAME,'Interp','Latex');
view(-37.5,30);

figure(2);
line([0 1],[0 0],'LineWidth',2,'Color','k');
line([0 0],[0 1],'LineWidth',3,'Color','k');
axis([0 1.3 0 1.3]);
axis off
axis square
text(1.06,0.01,0,'$x$','FontSize',FSIZE,'FontName',FNAME,'Interp','Latex');
text(-0.09,1.20,0,'$y$','FontSize',FSIZE,'FontName',FNAME,'Interp','Latex');

figure(3);
line([0 1],[0 0],'LineWidth',2,'Color','k');
line([1 1],[0 1],'LineWidth',2,'Color','k');
axis([0 1.3 0 1.3]);
axis off
axis square
text(-0.25,0.02,0,'$x$','FontSize',FSIZE,'FontName',FNAME,'Interp','Latex');
text(0.93,1.20,0,'$z$','FontSize',FSIZE,'FontName',FNAME,'Interp','Latex');

figure(4);
line([0 1],[0 0],'LineWidth',2,'Color','k');
line([0 0],[0 1],'LineWidth',3,'Color','k');
axis([0 1.3 0 1.3]);
axis off
axis square
text(1.05,-0.01,0,'$y$','FontSize',FSIZE,'FontName',FNAME,'Interp','Latex');
text(-0.07,1.18,0,'$z$','FontSize',FSIZE,'FontName',FNAME,'Interp','Latex');

% figure(1); print -dpng -r400 ../../figures/chapter01/fig_xyz_view1.png
% figure(2); print -dpng -r400 ../../figures/chapter01/fig_xyz_view2.png
% figure(3); print -dpng -r400 ../../figures/chapter01/fig_xyz_view3.png
% figure(4); print -dpng -r400 ../../figures/chapter01/fig_xyz_view4.png
