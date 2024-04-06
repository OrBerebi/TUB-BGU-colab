%FIG_ROTATION generates Figure 1.15, showing balloon plots of
% truncated spherical cap functions, illustrating rotation.
%
% Fundmentals of Spherical Array Processing
% Boaz Rafaely, 2017.

close all;
clear all;

path(path,'../../math');
path(path,'../../plot');

AxisFontSize=12;

% spherical cap
N=2;
alpha=(pi/180)*30;
fnm=zeros((N+1)^2,1);
fnm(1)=sqrt(pi)*(1-cos(alpha));
for n=1:N,
    fn0=sqrt(pi/(2*n+1))*(legendreP(n-1,cos(alpha))-legendreP(n+1,cos(alpha)));
    fnm(n^2+1:n^2+2*n+1)=[zeros(n,1);fn0;zeros(n,1)];
end;

% Rotations
fnm1=wignerD(N,0,pi/4,0)*fnm;
fnm2=wignerD(N,pi/4,0,0)*fnm1;
fnm3=wignerD(N,0,pi/4,0)*fnm2;


figure(1);

subplot(221);
plot_balloon(fnm,[180,0],0.05);
axis tight;
title('Original','Interp','Latex','FontSize',AxisFontSize);

subplot(222);
plot_balloon(fnm1,[180,0],0.05);
axis tight;
camzoom(0.85);
title('$\Lambda(0,45^\circ,0)$','Interp','Latex','FontSize',AxisFontSize);

subplot(223);
plot_balloon(fnm2,[180,0],0.05);
axis tight;
camzoom(0.85);
title('$\Lambda(45^\circ,45^\circ,0)$','Interp','Latex','FontSize',AxisFontSize);

subplot(224);
plot_balloon(fnm3,[180,0],0.05);
camzoom(0.9);
axis tight;
title('$\Lambda(45^\circ,90^\circ,0)$','Interp','Latex','FontSize',AxisFontSize);

% print -dpng -r400 ../../figures/chapter01/fig_rotation.png


