%FIG_RADIAL_FUNCTIONS_2 generates Figures 2,6 and 2.10, 
% illustrating the spherical Bessel function and the 
% rigid-sphere radial function as a function of the order n. 
%
% Fundmentals of Spherical Array Processing
% Boaz Rafaely, 2017.

close all;
clear all;

path(path,'../../math');
path(path,'../../plot');

AxisFontSize=14;

N=30;
kr16=16;
kr8=8;

for n=0:N,    
    j8(n+1)=Bn(n,kr8,kr8,0);
    j16(n+1)=Bn(n,kr16,kr16,0);
    b8(n+1)=Bn(n,kr8,kr8,1);
    b16(n+1)=Bn(n,kr16,kr16,1);    
end;

figure(1);
plot(0:N,20*log10(abs(j8)),'-o','LineWidth',2,'Color',[0 0 0.5]);
axis([0 N -100 20]);
xlabel('$n$','FontSize',AxisFontSize,'Interp','Latex');
ylabel('Magnitude (dB)','FontSize',AxisFontSize,'Interp','Latex');
title('$kr=8$','FontSize',AxisFontSize,'Interp','Latex');
set(gca,'FontSize',AxisFontSize);

figure(2);
plot(0:N,20*log10(abs(j16)),'-o','LineWidth',2,'Color',[0 0 0.5]);
axis([0 N -100 20]);
xlabel('$n$','FontSize',AxisFontSize,'Interp','Latex');
ylabel('Magnitude (dB)','FontSize',AxisFontSize,'Interp','Latex');
title('$kr=16$','FontSize',AxisFontSize,'Interp','Latex');
set(gca,'FontSize',AxisFontSize);

figure(3);
plot(0:N,20*log10(abs(b8)),'-o','LineWidth',2,'Color',[0 0 0.5]);
axis([0 N -100 20]);
xlabel('$n$','FontSize',AxisFontSize,'Interp','Latex');
ylabel('Magnitude (dB)','FontSize',AxisFontSize,'Interp','Latex');
title('$kr=8$','FontSize',AxisFontSize,'Interp','Latex');
set(gca,'FontSize',AxisFontSize);

figure(4);
plot(0:N,20*log10(abs(b16)),'-o','LineWidth',2,'Color',[0 0 0.5]);
axis([0 N -100 20]);
xlabel('$n$','FontSize',AxisFontSize,'Interp','Latex');
ylabel('Magnitude (dB)','FontSize',AxisFontSize,'Interp','Latex');
title('$kr=16$','FontSize',AxisFontSize,'Interp','Latex');
set(gca,'FontSize',AxisFontSize);

% figure(1); print -dpng -r400 ../../figures/chapter02/fig_rad_fun_n_jn1.png
% figure(2); print -dpng -r400 ../../figures/chapter02/fig_rad_fun_n_jn2.png
% figure(3); print -dpng -r400 ../../figures/chapter02/fig_rad_fun_n_bn1.png
% figure(4); print -dpng -r400 ../../figures/chapter02/fig_rad_fun_n_bn2.png
