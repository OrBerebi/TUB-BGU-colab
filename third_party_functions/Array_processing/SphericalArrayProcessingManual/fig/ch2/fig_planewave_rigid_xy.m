%FIG_PLANEWAVE_RIGID_XY generates figure 2.11, 
% illustrating a spherical-harmonics composition 
% of a plane wave around a rigid sphere, presented
% over the xy plane. 
%
% Fundmentals of Spherical Array Processing
% Boaz Rafaely, 2017.

close all;
clear all;

path(path,'../../math');
path(path,'../../plot');

AxisFontSize=14;

% Unit circle
z=linspace(0,2*pi,300);

N=32;
x=linspace(-20,20,100);
y=linspace(-20,20,100);

[X,Y]=meshgrid(x,y);
X1=reshape(X,length(x)*length(y),1);
Y1=reshape(Y,length(x)*length(y),1);
ph=atan2(Y1,X1);
th=pi/2*ones(size(ph));
r=sqrt(X1.^2+Y1.^2);

k=1;
thk=pi/2;
phk=pi/9;
ra1=10;
ra2=3;
ra3=1;

B11=BnMat(N,k*r',k*ra1,1);
B22=BnMat(N,k*r',k*ra2,1);
B33=BnMat(N,k*r',k*ra3,1);
Yk=conj(sh2(N,thk,phk));
Y=(sh2(N,th,ph)).';

p1=B11.*Y*Yk;
p2=B22.*Y*Yk;
p3=B33.*Y*Yk;

% zero region inside rigid sphere
for q=1:length(r),
    if r(q)<=ra1, p1(q)=0; end;
    if r(q)<=ra2, p2(q)=0; end;
    if r(q)<=ra3, p3(q)=0; end;
end;
        

p1a=reshape(p1,length(x),length(y));
p2a=reshape(p2,length(x),length(y));
p3a=reshape(p3,length(x),length(y));

figure(1);
[c,h]=contourf(x,y,real(p1a),'LineStyle','none'); 
axis square
colormap(jet);
colorbar;
set(gca,'FontSize',AxisFontSize);
set(colorbar,'FontSize',AxisFontSize);
xlabel('$x$ (m)','FontSize',AxisFontSize,'Interp','Latex');
ylabel('$y$ (m)','FontSize',AxisFontSize,'Interp','Latex');
title(strcat('$r_a=$',num2str(ra1)),'FontSize',AxisFontSize,'Interp','Latex');
hold on;
fill(ra1*cos(z),ra1*sin(z),'w');

figure(2);
[c,h]=contourf(x,y,real(p2a),'LineStyle','none'); 
axis square
colormap(jet);
colorbar;
set(gca,'FontSize',AxisFontSize);
xlabel('$x$ (m)','FontSize',AxisFontSize,'Interp','Latex');
ylabel('$y$ (m)','FontSize',AxisFontSize,'Interp','Latex');
title(strcat('$r_a=$',num2str(ra2)),'FontSize',AxisFontSize,'Interp','Latex');
hold on;
fill(ra2*cos(z),ra2*sin(z),'w');

figure(3);
[c,h]=contourf(x,y,real(p3a),'LineStyle','none'); 
axis square
colormap(jet);
colorbar;
set(gca,'FontSize',AxisFontSize);
xlabel('$x$ (m)','FontSize',AxisFontSize,'Interp','Latex');
ylabel('$y$ (m)','FontSize',AxisFontSize,'Interp','Latex');
title(strcat('$r_a=$',num2str(ra3)),'FontSize',AxisFontSize,'Interp','Latex');
hold on;
fill(ra3*cos(z),ra3*sin(z),'w');

% figure(1); print -dpng -r600 ../../figures/chapter02/fig_pressureNxy_rigid1.png
% figure(2); print -dpng -r600 ../../figures/chapter02/fig_pressureNxy_rigid2.png
% figure(3); print -dpng -r600 ../../figures/chapter02/fig_pressureNxy_rigid3.png
