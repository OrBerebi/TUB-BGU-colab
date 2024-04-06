%FIG_PN generates plots of the Legendre polynomial in Figures 1.11.
%
% Fundmentals of Spherical Array Processing
% Boaz Rafaely, 2017.

close all;
clear all;

AxisFontSize=14;

x=linspace(-1,1,256);
P0=legendreP(0,x);
P1=legendreP(1,x);
P2=legendreP(2,x);
P3=legendreP(3,x);
P4=legendreP(4,x);

figure(1);

subplot(3,2,1); 
plot(x,P0,'k-','LineWidth',1.5,'Color',[0 0 0.5]);
text(-0.95,2.3,'$n=0$','FontSize',AxisFontSize,'Interp','Latex'); 
set(gca,'FontSize',AxisFontSize,'XTickLabel','');

subplot(3,2,3); 
plot(x,P1,'k-','LineWidth',1.5,'Color',[0 0 0.5]); 
text(-0.95,1.3,'$n=1$','FontSize',AxisFontSize,'Interp','Latex'); 
set(gca,'FontSize',AxisFontSize,'XTickLabel','');

subplot(3,2,4); 
h=plot(x,P2,'k-','LineWidth',1.5,'Color',[0 0 0.5]); 
text(-0.95,1.3,'$n=2$','FontSize',AxisFontSize,'Interp','Latex');
set(gca,'FontSize',AxisFontSize,'XTickLabel','');

subplot(3,2,5); 
plot(x,P3,'k-','LineWidth',1.5,'Color',[0 0 0.5]); 
text(-0.95,1.3,'$n=3$','FontSize',AxisFontSize,'Interp','Latex'); 
set(gca,'FontSize',AxisFontSize);
xlabel('$x$','FontSize',AxisFontSize,'Interp','Latex');

subplot(3,2,6); 
plot(x,P4,'k-','LineWidth',1.5,'Color',[0 0 0.5]); 
text(-0.95,1.3,'$n=4$','FontSize',AxisFontSize,'Interp','Latex'); 
set(gca,'FontSize',AxisFontSize);
xlabel('$x$','FontSize',AxisFontSize,'Interp','Latex');

% print -dpng -r600 ../../figures/chapter01/fig_Pn.png


