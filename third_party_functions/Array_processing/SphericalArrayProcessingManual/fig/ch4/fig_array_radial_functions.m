%FIG_ARRAY_RADIAL_FUNCTIONS generates Figures 4.2, 4.5, 4.6, 
% showing radial functions to compare various array designs.
%
% Fundmentals of Spherical Array Processing
% Boaz Rafaely, 2017.

close all;
clear all;

path(path,'../../math');
path(path,'../../plot');

AxisFontSize=14;

N=3;
kr=linspace(0.01,6,10000);
ka=kr;
alpha=0.8333; 

BBo=zeros(N+1,length(kr)); % open
BBo2=zeros(N+1,length(kr)); % open 2
BBr=zeros(N+1,length(kr)); % rigid
BBc=zeros(N+1,length(kr)); % cardioid
BBd=zeros(N+1,length(kr)); % dual sphere

for n=0:N,    
    BBo(n+1,:)=Bn(n,kr,ka,0);
    BBo2(n+1,:)=Bn(n,kr*alpha,ka*alpha,0);
    BBr(n+1,:)=Bn(n,kr,ka,1);
    BBc(n+1,:)=Bn(n,kr,ka,2);
    
    beta = 0.5 + 0.5*sign( abs(Bn(n,kr*alpha,kr*alpha,0)) - abs(Bn(n,kr,kr,0)) );
    BBd(n+1,:) = (1-beta).*Bn(n,kr,kr,0) + beta.* Bn(n,kr*alpha,kr*alpha,0);
end;


figure;
for n=0:N,
    plot(kr,20*log10(abs(BBo(n+1,:))),'--','LineWidth',1.5,'Color',[0 0 0]); hold on;
    plot(kr,20*log10(abs(BBr(n+1,:))),'-','LineWidth',2,'Color',[0 0 0.5]); hold on;
end;
set(gca,'FontSize',18);
xlabel('$kr$','FontSize',AxisFontSize,'Interp','Latex');
ylabel('Magnitude (dB)','FontSize',AxisFontSize,'Interp','Latex');
axis([min(kr) max(kr) -60 30]);
legend('Open sphere','Rigid sphere');


figure;
for n=0:N,
    plot(kr,20*log10(abs(BBo(n+1,:))),'--','LineWidth',1.5,'Color',[0 0 0]); hold on;
    plot(kr,20*log10(abs(BBc(n+1,:))),'-','LineWidth',2,'Color',[0 0 0.5]); hold on;
end;    
hold off;
set(gca,'FontSize',18);
xlabel('$kr$','FontSize',AxisFontSize,'Interp','Latex');
ylabel('Magnitude (dB)','FontSize',AxisFontSize,'Interp','Latex');
axis([min(kr) max(kr) -60 30]);
legend('Pressure mic.','Cardioid mic.');


figure;
for n=0:N,
    plot(kr,20*log10(abs(BBo(n+1,:))),'-','LineWidth',1,'Color',[0 0 0]); hold on;
    plot(kr,20*log10(abs(BBo2(n+1,:))),'--','LineWidth',1.5,'Color',[0 0 0.5]); hold on;
    plot(kr,20*log10(abs(BBd(n+1,:))),'-','LineWidth',2.5,'Color',[0 0 0.5]); hold on;
end;
hold off;
legend('Open r_1','Open r_2','Dual sphere','Location','SouthWest');
set(gca,'FontSize',18);
xlabel('$k$','FontSize',AxisFontSize,'Interp','Latex');
ylabel('Magnitude (dB)','FontSize',AxisFontSize,'Interp','Latex');
axis([min(kr) max(kr) -60 30]);

% figure(1); print -dpng -r400 ../../figures/chapter04/fig_bn_rigid.png
% figure(2); print -dpng -r400 ../../figures/chapter04/fig_bn_card.png
% figure(3); print -dpng -r400 ../../figures/chapter04/fig_bn_dual.png
