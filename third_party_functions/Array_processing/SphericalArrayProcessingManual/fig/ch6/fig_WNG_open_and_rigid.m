% FIG_WNG_OPEN_AND_RIGID generates Figure 6.2, 
% illustrating white-noise gain for
% open and rigid sphere arrays.
%
% Fundmentals of Spherical Array Processing
% Boaz Rafaely, 2017.

close all;
clear all;

path(path,'../../math');
path(path,'../../plot');

AxisFontSize=14;

N=3;
[a,th,ph]=uniform_sampling(N);
Q=length(a);
kr=linspace(0,N,512);
Y0=sh2(N,0,0);

Bopen=BnMat(N,kr,kr,0);
vnm=Bopen*diag(conj(Y0));
WNGopen=(Q/(4*pi))*sum(abs(vnm).^2,2);

Brigid=BnMat(N,kr,kr,1);
vnm=Brigid*diag(conj(Y0));
WNGrigid=(Q/(4*pi))*sum(abs(vnm).^2,2);

figure;
plot(kr,10*log10(WNGopen),'-','LineWidth',2,'Color',[0 0 0.5]); hold on;
plot(kr,10*log10(WNGrigid),'-','LineWidth',1,'Color',[0 0 0]);
plot(kr,ones(size(kr))*10*log10(Q),'k--','LineWidth',2,'Color',[0 0 0.5]);
lg=legend('Open','Rigid','Q');
set(lg,'FontSize',AxisFontSize);
axis([0,3,14,19]);
xlabel('$kr$','FontSize',AxisFontSize,'Interp','Latex');
ylabel('$WNG\,\,$ (dB)','FontSize',AxisFontSize,'Interp','Latex');
set(gca,'FontSize',AxisFontSize);

% print -dpng -r400 ../../figures/chapter06/fig_WNG.png

