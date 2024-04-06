%FIG_SINC_AND_CAP generates Figures 1.12 and 1.13, illustrating a truncated
% spherical harmonics series function and a spherical cap function. 
%
% Fundmentals of Spherical Array Processing
% Boaz Rafaely, 2017.

close all;
clear all;

path(path,'../../math');
path(path,'../../plot');

AxisFontSize=14;

% Truncated shperical harmonics series
N=[8 20];
Theta=linspace(-pi/2,pi/2,500);
z=cos(Theta);
S1=zeros(length(N),length(z));

for i=1:length(N),
    n=N(i);
    S1(i,:)=((n+1)./(4*pi*(z-1))) .* (legendreP(n+1,z)-legendreP(n,z)) ;
end;

% Spherical cap
N=20;
alpha=(pi/180)*[15 45];
S2=zeros(length(alpha),N+1);
S2(:,1)=sqrt(pi)*(1-cos(alpha)).';

for n=1:N,
    S2(:,n+1)=sqrt(pi/(2*n+1)) * (legendreP(n-1,cos(alpha))-legendreP(n+1,cos(alpha))).' ;
end;

TH=linspace(0,pi/2,200);
f1=(TH<=alpha(1));
f2=(TH<=alpha(2));


figure(1);
plot(Theta*180/pi,S1(1,:),'--','LineWidth',1.5,'Color',[0 0 0]);
hold on;
plot(Theta*180/pi,S1(2,:),'-','LineWidth',2,'Color',[0 0 0.5]);
set(gca,'FontSize',AxisFontSize);
ylabel('Amplitude','Interp','Latex','FontSize',AxisFontSize);
xlabel('$\Theta$ (degrees)','Interp','Latex','FontSize',AxisFontSize);
hold off;
lg=legend('N=8','N=20');
set(lg,'FontSize',AxisFontSize);
axis([-90,90,-10,40]);

figure(2);
plot(0:N,S2(1,:),'x--','LineWidth',1.5,'Color',[0 0 0],'MarkerSize',10);
hold on;
plot(0:N,S2(2,:),'o-','LineWidth',2,'Color',[0 0 0.5],'MarkerSize',8);
set(gca,'FontSize',AxisFontSize);
ylabel('Amplitude','Interp','Latex','FontSize',AxisFontSize);
xlabel('$n$','Interp','Latex','FontSize',AxisFontSize);
hold off;
lg=legend('\alpha=15^\circ','\alpha=45^\circ');
set(lg,'FontSize',AxisFontSize);

figure(3);
plot(TH*180/pi,f1,'--','LineWidth',1.5,'Color',[0 0 0]);
hold on;
plot(TH*180/pi,f2,'-','LineWidth',2,'Color',[0 0 0.5]);
set(gca,'FontSize',AxisFontSize);
ylabel('Amplitude','Interp','Latex','FontSize',AxisFontSize);
xlabel('$\theta$ (degrees)','Interp','Latex','FontSize',AxisFontSize);
hold off;
lg=legend('\alpha=15^\circ','\alpha=45^\circ');
set(lg,'FontSize',AxisFontSize);
axis([0,90,0,1.1]);

% figure(1); print -dpng -r400 ../../figures/chapter01/fig_sincSH1.png
% figure(2); print -dpng -r400 ../../figures/chapter01/fig_sincSH2.png
% figure(3); print -dpng -r400 ../../figures/chapter01/fig_sincSH3.png
