%FIG_DOLPHCHEBYSHEV_BEAMPATTERN generates Figures 6.7, 6.8, 
% illustrating Dolph-Chebyshev beam patterns.
%
% Fundmentals of Spherical Array Processing
% Boaz Rafaely, 2017.

clear all
close all;
clc;

path(path,'../../math');
path(path,'../../plot');

AxisFontSize=14;

M=4001; % points for plot 
theta = linspace(-pi,pi,M);


% First plot

R=10^(25/20); % side lobe level constraint in dB

N=4; 
x0=cosh(1/(2*N)*acosh(R));
theta0=2*acos((1/x0)*cos(pi/(4*N)))*180/pi;
dn=(2*pi/R)*DolphPACT(N)*(x0.^(2*(0:N))).';
y1a=B(theta,N,dn);

N=9; 
x0=cosh(1/(2*N)*acosh(R));
theta0=2*acos((1/x0)*cos(pi/(4*N)))*180/pi;
dn=(2*pi/R)*DolphPACT(N)*(x0.^(2*(0:N))).';
y1b=B(theta,N,dn);

figure;
plot(theta*180/pi,20*log10(abs(y1a)),'-','LineWidth',2,'Color',[0 0 0.5]); hold on;
plot(theta*180/pi,20*log10(abs(y1b)),'--','LineWidth',1.5,'Color',[0 0 0]);
set(gca,'FontSize',AxisFontSize);
xlabel('$\theta$ (degrees)','interp','Latex','FontSize',AxisFontSize);
ylabel('Magnitude (dB)','interp','Latex','FontSize',AxisFontSize);
lg=legend('N=4','N=9');
set(lg,'FontSize',AxisFontSize);
axis([-180 180 -35 5]);



% second plot

theta0=45*pi/180;

N=4; 
x0=cos(pi/(4*N))/(cos(theta0/2));
R=cosh(2*N*acosh(x0));
dn=(2*pi/R)*DolphPACT(N)*(x0.^(2*(0:N))).';
y2a=B(theta,N,dn);

N=9; 
x0=cos(pi/(4*N))/(cos(theta0/2));
R=cosh(2*N*acosh(x0));
dn=(2*pi/R)*DolphPACT(N)*(x0.^(2*(0:N))).';
y2b=B(theta,N,dn);

figure;
plot(theta*180/pi,20*log10(abs(y2a)),'-','LineWidth',2,'Color',[0 0 0.5]); hold on;
plot(theta*180/pi,20*log10(abs(y2b)),'--','LineWidth',1.5,'Color',[0 0 0]);
set(gca,'FontSize',AxisFontSize);
xlabel('$\theta$ (degrees)','interp','Latex','FontSize',AxisFontSize)
ylabel('Magnitude (dB)','interp','Latex','FontSize',AxisFontSize);
lg=legend('N=4','N=9');
set(lg,'FontSize',AxisFontSize);
axis([-180 180 -80 5]);

% figure(1); print -dpng -r400 ../../figures/chapter06/fig_dolph_cheby1.png
% figure(2); print -dpng -r400 ../../figures/chapter06/fig_dolph_cheby2.png



function Y = B(theta,N,dn)
% calculates the transform of a window

for n=0:N
    L(n+1,:) = dn(n+1)*((2*n+1)/(4*pi)) * legendreP(n,cos(theta));
end
Y=sum(L,1);
end
