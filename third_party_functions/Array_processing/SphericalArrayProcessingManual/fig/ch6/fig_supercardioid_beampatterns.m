%FIG_SUPERCARDIOID_BEAMPATTERNS generates Figure 6.4, 
% illustrating super-cardioid beam patterns.
%
% Fundmentals of Spherical Array Processing
% Boaz Rafaely, 2017.

close all;
clear all;

path(path,'../../math');
path(path,'../../plot');

AxisFontSize=14;

% legendre polynomials
p0=[1 0 0 0 0];
p1=[0 1 0 0 0 ];
p2=0.5*[-1 0 3 0 0];
p3=0.5*[0 -3 0 5 0];
p4=(1/8)*[3 0 -30 0 35];
p=[p0;p1;p2;p3;p4];


for N=0:4,

A=0;
B=0;
for n=0:N,
    for nn=0:N,
        a=0;
        b=0;
        for k=0:n,
            for l=0:nn,
                a=a+(1/(k+l+1))*p(n+1,k+1)*p(nn+1,l+1);
                b=b+(((-1)^(k+l))/(k+l+1))*p(n+1,k+1)*p(nn+1,l+1);
            end;
        end;
        A(n+1,nn+1)=(1/(8*pi))*(2*n+1)*(2*nn+1)*a;
        B(n+1,nn+1)=(1/(8*pi))*(2*n+1)*(2*nn+1)*b;
    end;
end;

[V,D]=eig(A,B,'chol');
[lambda,I]=max(diag(D));
dn=V(:,I);

pp=diag(dn)*diag((2*[0:N]+1)/(4*pi))*p(1:N+1,:);
yp=sum(pp);

Theta=linspace(0,2*pi,512);
x=cos(Theta);
y=0;
for i=0:length(yp)-1,
    y=y+yp(i+1)*x.^i;
end;
y=y/abs(y(1))-1e-3;



F(N+1)=10*log10(lambda);
Y(N+1,:)=y;

end;



figure;

subplot(221);
h1=polarplot(Theta,abs(Y(2,:)),'-');
set(h1,'LineWidth',2,'Color',[0 0 0.5]);
set(gca,'FontSize',AxisFontSize-2);
title('$N=1, F=11.4\,$dB','FontSize',AxisFontSize,'interp','Latex');

subplot(222);
h1=polarplot(Theta,abs(Y(3,:)),'-');
set(h1,'LineWidth',2,'Color',[0 0 0.5]);
set(gca,'FontSize',AxisFontSize-2);
title('$N=2, F=24.1\,$dB','FontSize',AxisFontSize,'interp','Latex');

subplot(223);
h1=polarplot(Theta,abs(Y(4,:)),'-');
set(h1,'LineWidth',2,'Color',[0 0 0.5]);
set(gca,'FontSize',AxisFontSize-2);
title('$N=3, F=37.7\,$dB','FontSize',AxisFontSize,'interp','Latex');

subplot(224);
h1=polarplot(Theta,abs(Y(5,:)),'-');
set(h1,'LineWidth',2,'Color',[0 0 0.5]);
set(gca,'FontSize',AxisFontSize-2);
title('$N=4, F=51.8\,$dB','FontSize',AxisFontSize,'interp','Latex');

% print -dpng -r600 ../../figures/chapter06/fig_super_cardioid.png

