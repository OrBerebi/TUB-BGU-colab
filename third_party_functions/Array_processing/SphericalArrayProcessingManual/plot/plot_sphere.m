function plot_sphere(fnm);
%PLOT_SPHERE generates a plot of a function on a sphere.
% plot_sphere(fnm);
% fnm is the spherical harmonics coefficient vector.
% The function plots real{f}.
%
% Fundmentals of Spherical Array Processing
% Boaz Rafaely, 2017.

N=sqrt(length(fnm))-1;

AxisFontSize=14;

n=150;
[x,y,z] = sphere(n);
[TH,PH,R]=c2s(x,y,z);
th=reshape(TH,size(TH,1)*size(TH,2),1);
ph=reshape(PH,size(PH,1)*size(PH,2),1);

Y=(sh2(N,th,ph)).';
f=Y*fnm;
F=reshape(f,size(TH,1),size(TH,2));
F=real(F);

h1=surf(x,y,z,F);

colormap(jet);
brighten(0.3);
set(h1,'EdgeAlpha',0.01);

colorbar;
caxis([0.1*floor(min(min(10*F))),0.1*ceil(max(max(10*F)))]);
set(colorbar,'FontSize',AxisFontSize);

axis tight equal
xlabel('$x$','FontSize',AxisFontSize,'Interp','Latex');
ylabel('$y$','FontSize',AxisFontSize,'Interp','Latex');
zlabel('$z$','FontSize',AxisFontSize,'Interp','Latex');
set(gca,'FontSize',AxisFontSize);


