%FIG_MULTIPLE_OBJECTIVE_BEAMPATTERNS generates Figures 6.9, 6.10, 
% illustrating beam patterns for array design with multiple 
% objectives and constraints.
% It uses the Matlab optimization function fmincon.
%
% Fundmentals of Spherical Array Processing
% Boaz Rafaely, 2017.

close all;
clear all;
clc;

path(path,'../../math');
path(path,'../../plot');

AxisFontSize=14;

N=4;
kr=2;
sphere=1; % rigid
[a,th,ph]=uniform_sampling(N);
Q=length(a);
nn=(2*[0:N]+1).';
vn=(1/(4*pi))*nn;
bn=BnMat(N,kr,kr,sphere); 
bn=bn([1:N+1].^2).';
A = (4*pi/Q)*diag(vn)*diag(1./abs(bn).^2);
B = (1/(4*pi))*diag(vn);
WNGmin=10;
l2SL=10^(-30/10);

THi=(pi/180)*linspace(60,180,50);
for n=0:N,
    P=legendre(n,cos(THi));
    Pn2(n+1,:)=P(1,:);
end;
Vn=diag(vn)*Pn2;

x0=zeros(N+1,1);

% design 1, robust max directivity
dn1 = fmincon(@(x)(x'*B*x),x0,[],[],vn',1,[],[],@(x)deal(x'*A*x-1/WNGmin,0));

% design 2, robust max side-lobe level
dn2 = fmincon(@(x)(x'*B*x),x0,[],[],vn',1,[],[],@(x)deal([x'*A*x-1/WNGmin,abs(x'*Vn).^2-l2SL],0));


DI1=abs(dn1'*vn)^2/(dn1'*B*dn1);
WNG1=abs(dn1'*vn)^2/(dn1'*A*dn1);
SL1=max(abs(dn1'*Vn));

DI2=abs(dn2'*vn)^2/(dn2'*B*dn2);
WNG2=abs(dn2'*vn)^2/(dn2'*A*dn2);
SL2=max(abs(dn2'*Vn));

% plot
thp=(pi/180)*linspace(-180,180,1024);
for n=0:N,
    P=legendre(n,cos(thp));
    Pn(n+1,:)=P(1,:);
end;
y1=dn1.'*diag(vn)*Pn;
y2=dn2.'*diag(vn)*Pn;

figure;
plot((180/pi)*thp,20*log10(abs(y1)),'-','LineWidth',2,'Color',[0 0 0.5]);
xlabel('$\Theta$ (degrees)','FontSize',AxisFontSize,'interp','Latex');
ylabel('Magnitude (dB)','FontSize',AxisFontSize,'interp','Latex');
set(gca,'FontSize',AxisFontSize);
axis([-180,180,-70,10]);

figure;
plot((180/pi)*thp,20*log10(abs(y2)),'k-','LineWidth',2,'Color',[0 0 0.5]);
xlabel('$\Theta$ (degrees)','FontSize',AxisFontSize,'interp','Latex');
ylabel('Magnitude (dB)','FontSize',AxisFontSize,'interp','Latex');
set(gca,'FontSize',AxisFontSize);
axis([-180,180,-70,10]);

fprintf('\nDF1=%1.2f',DI1);
fprintf('\nWNG1=%1.2f',WNG1);
fprintf('\nSL1=%1.2f',20*log10(SL1));
fprintf('\nDF2=%1.2f',DI2);
fprintf('\nWNG2=%1.2f',WNG2);
fprintf('\nSL2=%1.2f\n',20*log10(SL2));

% figure(1); print -dpng -r400 ../../figures/chapter06/fig_multiple1.png
% figure(2); print -dpng -r400 ../../figures/chapter06/fig_multiple2.png


