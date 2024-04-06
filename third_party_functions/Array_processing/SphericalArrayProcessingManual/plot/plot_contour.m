function y=plot_contour(fnm,normalization,absolute);
%PLOT_CONTOUR generates a contour plot.
% y = plot_contour(fnm,normalization,absolute);
% fnm is the spherical harmonics coefficient vector.
% normalization (0/1) is a normalisation option (deafault=0).
% absolute (0/1) applies absolute value before plotting (deafault=1).
%
% Fundmentals of Spherical Array Processing
% Boaz Rafaely, 2017.

if nargin<3 | isempty(absolute), absolute=1; end; % deafault - apply abs() to function
if nargin<2 | isempty(normalization), normalization=0; end; % deafault - no amplitude normalization 

N=sqrt(length(fnm))-1;

AxisFontSize=14;

Np=29;
[ap,thp,php]=equiangle_sampling(Np);
Yp=sh2(N,thp,php); 
f=fnm.'*Yp;

% Amplitude normalization
if normalization
    f=f/max(abs(f));
end

% Absolute value
if absolute
    f=abs(f);
end

% grid
phi=php(1:2*(Np+1));
the=thp(1:2*(Np+1):end);

% reshape to a matrix
F=(reshape(f,2*(Np+1),2*(Np+1)).');

% ignore small imaginary part
if max(max(imag(F)))<1e-10, F=real(F); end

% plot
[c,h]=contourf(phi*180/pi,the*180/pi,F); 
colormap(flipud(bone));
colorbar;
caxis([0.1*floor(min(min(10*F))),0.1*ceil(max(max(10*F)))]);
set(colorbar,'FontSize',AxisFontSize);

set(gca,'FontSize',AxisFontSize);
xlabel('$\phi\,$ (degrees)','FontSize',AxisFontSize,'Interp','Latex');
ylabel('$\theta\,$ (degrees)','FontSize',AxisFontSize,'Interp','Latex');
