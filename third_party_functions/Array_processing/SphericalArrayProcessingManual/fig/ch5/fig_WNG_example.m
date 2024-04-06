%FIG_WNG_EXAMPLE generates Figure 5.5, 
% illustrating array white-noise gain.
%
% Fundmentals of Spherical Array Processing
% Boaz Rafaely, 2017.

close all;
clear all;

path(path,'../../math');
path(path,'../../plot');

AxisFontSize=14;

Q=9; N=2;
n=(0:N).';
vn=(1/(4*pi))*(2*n+1);
kr=linspace(0.01,N,512);
BB=BnMat(N,kr,kr,0);
Bn=BB(:,(n+1).^2);
dn=ones(N+1,1);
A=vn*vn.';

for i=1:length(kr),
    B=(4*pi/Q)*diag(vn)*diag(1./abs(Bn(i,:)).^2);
    WNG(i)=(dn'*A*dn)/(dn'*B*dn);
end;

h=plot(kr,WNG,'-','LineWidth',2,'Color',[0 0 0.5]); hold on
plot(kr,ones(size(kr)),'--','LineWidth',1.5,'Color',[0 0 0])
xlabel('$kr$','FontSize',AxisFontSize,'Interp','Latex');
ylabel('$WNG$','FontSize',AxisFontSize,'Interp','Latex');
legend('Q=9,N=2','Q=1,N=0');
set(gca,'FontSize',AxisFontSize);

% print -dpng -r400 ../../figures/chapter05/fig_WNG_example.png

