%FIG_ARRAY_DESIGN_EXAMPLES generates Figures 4.1, 4.3, 
% showing open-sphere and rigid-sphere radial functions
% to illustrate array designs.
%
% Fundmentals of Spherical Array Processing
% Boaz Rafaely, 2017.

close all;
clear all;

path(path,'../../math');
path(path,'../../plot');

AxisFontSize=14;

N=8;
N0=5;
r=0.08;
kr=linspace(0.01,7,10000);
c=343;
kr2f=c/(2*pi*r);
f=kr*kr2f;

bn0=zeros(N+1,length(kr)); % open sphere
bn1=zeros(N+1,length(kr)); % rigid sphere
for n=0:N,    
    bn0(n+1,:)=Bn(n,kr,kr,0);
    bn1(n+1,:)=Bn(n,kr,kr,1);
end;


figure;
plot(f,20*log10(abs(bn0)),'-','LineWidth',1.5,'Color',[0 0 0.5]);
text(0.2*kr2f,25,'0','FontSize',AxisFontSize);
text(0.6*kr2f,12.2,'1','FontSize',AxisFontSize);
text(0.9*kr2f,3.1,'2','FontSize',AxisFontSize);
text(1.2*kr2f,-7.0,'3','FontSize',AxisFontSize);
text(1.6*kr2f,-15.7,'4','FontSize',AxisFontSize);
text(1.9*kr2f,-24.5,'5','FontSize',AxisFontSize);
text(2.2*kr2f,-33,'6','FontSize',AxisFontSize);
text(2.6*kr2f,-42,'7','FontSize',AxisFontSize);
text(2.9*kr2f,-50,'8','FontSize',AxisFontSize);
set(gca,'FontSize',AxisFontSize);
xlabel('Frequency (Hz)','FontSize',AxisFontSize,'Interp','Latex');
ylabel('Magnitude (dB)','FontSize',AxisFontSize,'Interp','Latex');
axis([min(f) max(f) -60 30]);
h=line([N0*kr2f,N0*kr2f],[-60,30]);
set(h,'Color','k','LineStyle','--','LineWidth',1.5');
text(1.03*N0*kr2f,20,'$kr=5$','FontSize',AxisFontSize,'Color','k','Interp','Latex');
hold on; 
g=plot(f(2083),20*log10(abs(bn0(6,2083))),'kx');
set(g,'LineWidth',2,'MarkerSize',AxisFontSize+2);

figure;
plot(f,20*log10(abs(bn1)),'-','LineWidth',1.5,'Color',[0 0 0.5]);
text(0.2*kr2f,25,'0','FontSize',AxisFontSize);
text(0.6*kr2f,12.2+4,'1','FontSize',AxisFontSize);
text(0.9*kr2f,3.1+4,'2','FontSize',AxisFontSize);
text(1.2*kr2f,-7.0+4,'3','FontSize',AxisFontSize);
text(1.6*kr2f,-15.7+4,'4','FontSize',AxisFontSize);
text(1.9*kr2f,-24.5+4,'5','FontSize',AxisFontSize);
text(2.2*kr2f,-33+4,'6','FontSize',AxisFontSize);
text(2.6*kr2f,-42+4,'7','FontSize',AxisFontSize);
text(2.9*kr2f,-50+4,'8','FontSize',AxisFontSize);
set(gca,'FontSize',AxisFontSize);
xlabel('Frequency (Hz)','FontSize',AxisFontSize,'Interp','Latex');
ylabel('Magnitude (dB)','FontSize',AxisFontSize,'Interp','Latex');
axis([min(f) max(f) -60 30]);
h=line([N0*kr2f,N0*kr2f],[-60,30]);
set(h,'Color','k','LineStyle','--','LineWidth',1.5');
text(1.03*N0*kr2f,20,'$kr=5$','FontSize',AxisFontSize,'Color','k','Interp','Latex');
hold on; 
g=plot(f(2083),20*log10(abs(bn1(6,2083))),'kx');
set(g,'LineWidth',2,'MarkerSize',AxisFontSize+2);

% figure(1); print -dpng -r400 ../../figures/chapter04/fig_design_examples_open.png
% figure(2); print -dpng -r400 ../../figures/chapter04/fig_design_examples_rigid.png
