%FIG_HYPERCARDIOID_BEAMPATTERN generates Figure 6.1, 
% illustrating hyper-cardioid beam patterns
%
% Fundmentals of Spherical Array Processing
% Boaz Rafaely, 2017.

close all;
clear all;

path(path,'../../math');
path(path,'../../plot');

AxisFontSize=14;

Theta=linspace(0,2*pi,512);
z=cos(Theta);

y0=ones(size(Theta));
y1=(1/4)*(3*z+1);
y2=(1/6)*(5*z.^2+2*z-1);
y3=(1/32)*(35*z.^3+15*z.^2-15*z-3);
y4=(1/40)*(63*z.^4+28*z.^3-42*z.^2-12*z+3);
y5=(1/96)*(231*z.^5+105*z.^4-210*z.^3-70*z.^2+35*z+5);

figure;
subplot(221);
h1=polarplot(Theta,(abs(y1)),'-');
set(h1,'LineWidth',2,'Color',[0 0 0.5]);
set(gca,'FontSize',AxisFontSize-2);
title('$N=1$','FontSize',AxisFontSize,'interp','Latex');

subplot(222);
h1=polarplot(Theta,(abs(y2)),'-');
set(h1,'LineWidth',2,'Color',[0 0 0.5]);
set(gca,'FontSize',AxisFontSize-2);
title('$N=2$','FontSize',AxisFontSize,'interp','Latex');

subplot(223);
h1=polarplot(Theta,(abs(y3)),'-');
set(h1,'LineWidth',2,'Color',[0 0 0.5]);
set(gca,'FontSize',AxisFontSize-2);
title('$N=3$','FontSize',AxisFontSize,'interp','Latex');

subplot(224);
h1=polarplot(Theta,(abs(y4)),'-');
set(h1,'LineWidth',2,'Color',[0 0 0.5]);
set(gca,'FontSize',AxisFontSize-2);
title('$N=4$','FontSize',AxisFontSize,'interp','Latex');

% print -dpng -r600 ../../figures/chapter06/fig_hyper_cardioid.png

