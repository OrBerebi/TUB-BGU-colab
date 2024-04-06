function y = plot_balloon(fnm,viewangle,transparency);
%PLOT_BALLOON generates a balloon plot.
% y = plot_balloon(fnm,viewangle,transparency)
% fnm is the spherical harmonics coefficient vector.
% viewang controls the view angle with azimuth=viewang(1) and
% elevation=viewang(2).
% transparency is the grid transparency level.
%
% Fundmentals of Spherical Array Processing
% Boaz Rafaely, 2017.

if nargin<3 | isempty(transparency), transparency=0.01; end;
if nargin<2 | isempty(viewangle), viewangle=[-37.5,30]; end;

N=sqrt(length(fnm))-1;

n=100;
[x,y,z] = sphere(n);
[th,ph,r]=c2s(x,y,z);
thp=reshape(th,size(th,1)*size(th,2),1);
php=reshape(ph,size(th,1)*size(th,2),1);
Yp=sh2(N,thp,php);

f=fnm.'*Yp;
F=reshape(f,size(th,1),size(th,2));
[X,Y,Z]=s2c(th,ph,abs(F));

h=surf(X,Y,Z,sign(real(F)));
axis equal
axis off;

colormap(flipud(cool));
brighten(0.9);
shading faceted;
set(h,'EdgeAlpha',transparency);
caxis([-1,1]);

view(viewangle);

lighting phong
camlight('right');



