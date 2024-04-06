%FIG_WNG_AND_DI_EXAMPLE generates Figure 6.3, 
% comparing direcitvity index and white-noise gain 
% of two designs.
%
% Fundmentals of Spherical Array Processing
% Boaz Rafaely, 2017.

close all;
clear all;

path(path,'../../math');
path(path,'../../plot');

AxisFontSize=14;

N=4;
[a,th,ph]=uniform_sampling(N);
Q=length(a);
kr=linspace(0,N+1,1024);
Y0=sh2(N,0,0);

bn=zeros(N+1,length(kr));
for n=0:N,
    bn(n+1,:)=4*pi*(i^n)*besseljs(n,kr);
end;

dn_maxDI=ones(N+1,length(kr));
DI_maxDI=abs(sum(diag(2*[0:N]+1)*dn_maxDI,1).^2)./sum(diag(2*[0:N]+1)*abs(dn_maxDI.^2),1);
WNG_maxDI=(Q/((4*pi)^2))*abs(sum(diag(2*[0:N]+1)*dn_maxDI,1).^2)./sum(diag(2*[0:N]+1)*abs((dn_maxDI./bn).^2),1);

dn_maxWNG=abs(bn).^2;
DI_maxWNG=abs(sum(diag(2*[0:N]+1)*dn_maxWNG,1).^2)./sum(diag(2*[0:N]+1)*abs(dn_maxWNG.^2),1);
WNG_maxWNG=(Q/((4*pi)^2))*abs(sum(diag(2*[0:N]+1)*dn_maxWNG,1).^2)./sum(diag(2*[0:N]+1)*abs((dn_maxWNG./bn).^2),1);

figure;
plot(kr,10*log10(DI_maxDI),'-','LineWidth',2,'Color',[0 0 0.5]); hold on;
plot(kr,10*log10(DI_maxWNG),'--','LineWidth',1.5,'Color',[0 0 0]);
lg=legend('Max DI','Max WNG');
set(lg,'FontSize',AxisFontSize);
axis([0,max(kr),0,20]);
xlabel('$kr$','FontSize',AxisFontSize,'Interp','Latex');
ylabel('$DI\,$ (dB)','FontSize',AxisFontSize,'Interp','Latex');
set(gca,'FontSize',AxisFontSize);

figure;
plot(kr,10*log10(WNG_maxDI),'-','LineWidth',2,'Color',[0 0 0.5]); hold on;
plot(kr,10*log10(WNG_maxWNG),'--','LineWidth',1.5,'Color',[0 0 0]);
lg=legend('Max DI','Max WNG');
set(lg,'FontSize',AxisFontSize);
axis([0,max(kr),-30,30]);
xlabel('$kr$','FontSize',AxisFontSize,'Interp','Latex');
ylabel('$WNG\,\,$ (dB)','FontSize',AxisFontSize,'Interp','Latex');
set(gca,'FontSize',AxisFontSize);

% figure(1); print -dpng -r400 ../../figures/chapter06/fig_design_example1.png
% figure(2); print -dpng -r400 ../../figures/chapter06/fig_design_example2.png
