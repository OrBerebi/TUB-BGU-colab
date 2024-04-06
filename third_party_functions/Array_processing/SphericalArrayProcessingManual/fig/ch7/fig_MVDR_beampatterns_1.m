%FIG_MVDR_BEAMPATTERNS_1 generates Figures 7.1, 7.2, 
% showing beam patterns of MVDR beamformers.
%
% Fundmentals of Spherical Array Processing
% Boaz Rafaely, 2017.

close all;
clear all;

path(path,'../../math');
path(path,'../../plot');

AxisFontSize=14;

% array parameter setting
N=4;
Q=36; 
kr=N;
bn=BnMat(N,kr,kr,1); % rigid sphere
B=diag(bn);
th0=pi/3; ph0=2*pi/10; % look direction
Y0=sh2(N,th0,ph0);
vnm0=B*conj(Y0);


% MVDR with sensor noise
wnm=vnm0/(vnm0'*vnm0);
wnmB=conj(B)*wnm; 
% Note: this leads to wnmB*Ynm=conj(y).
% Then plot abs(conj(y))=abs(y).
% This way plot_contour and plot_balloon can be used directly
                  
figure;                  
plot_contour(wnmB);
hold on; 
plot((180/pi)*ph0,(180/pi)*th0,'w+','MarkerSize',14,'LineWidth',3.0);

figure;
plot_balloon(wnmB,[40,30],0.1);
axis on; axis tight; axis equal;
xlabel('$x$','FontSize',AxisFontSize,'Interp','Latex');
ylabel('$y$','FontSize',AxisFontSize,'Interp','Latex');
zlabel('$z$','FontSize',AxisFontSize,'Interp','Latex');
set(gca,'FontSize',AxisFontSize);


% MVDR with one disturbance and sensor noise
sigma2s=1; % desired signal 
sigma2d=0.5; % disturbance
sigma2n=0.1; % sensor noise
th1=pi/3; ph1=2*pi*8/9; % disturbance direction
Y1=sh2(N,th1,ph1);
vnm1=B*conj(Y1);
Sxx = sigma2s*vnm0*vnm0' + sigma2d*vnm1*vnm1' + sigma2n*(4*pi/Q)*eye((N+1)^2);
wnm = inv(Sxx)'*vnm0/(vnm0'*inv(Sxx)*vnm0);
wnmB=conj(B)*wnm; 
% Note: this computes conj(y) and not y, see above

figure;
plot_contour(wnmB);
hold on; 
plot((180/pi)*ph0,(180/pi)*th0,'w+','MarkerSize',14,'LineWidth',3.0);
plot((180/pi)*ph1,(180/pi)*th1,'k+','MarkerSize',14,'LineWidth',3.0);


figure;
plot_balloon(wnmB,[40,30],0.1); 
axis on; axis tight; axis equal;
xlabel('$x$','FontSize',AxisFontSize,'Interp','Latex');
ylabel('$y$','FontSize',AxisFontSize,'Interp','Latex');
zlabel('$z$','FontSize',AxisFontSize,'Interp','Latex');
set(gca,'FontSize',AxisFontSize);

% figure(1); print -dpng -r400 ../../figures/chapter07/fig_mvdr1_contour1.png
% figure(2); print -dpng -r600 ../../figures/chapter07/fig_mvdr1_balloon1.png
% figure(3); print -dpng -r400 ../../figures/chapter07/fig_mvdr1_contour2.png
% figure(4); print -dpng -r600 ../../figures/chapter07/fig_mvdr1_balloon2.png

