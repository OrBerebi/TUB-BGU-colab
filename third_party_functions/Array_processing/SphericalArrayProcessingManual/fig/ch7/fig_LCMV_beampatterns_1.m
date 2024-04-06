%FIG_LCMV_BEAMPATTERNS_1 generates Figures 7.5-7.7, 
% showing beam patterns of LCMV beamformers. 
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
Q=36; % no. of samples in uniform sampling
kr=N;
bn=BnMat(N,kr,kr,1); % rigid sphere
B=diag(bn);

% desired signal
th0=60*pi/180; ph0=36*pi/180; % look direction
Y0=sh2(N,th0,ph0);
sigma2s=1; % desired signal 
vnm0=B*conj(Y0);

% disturbance / null
th1=60*pi/180; ph1=320*pi/180; % disturbance direction
Y1=sh2(N,th1,ph1);
vnm1=B*conj(Y1);

% Sxx and sensor noise
sigma2n=0.1;
Sxx = sigma2s*vnm0*vnm0' + sigma2n*(4*pi/Q)*eye((N+1)^2);

% constraints
Vnm=[vnm0,vnm1];
c=[1;0]


% LCMV with one null
wnm_conj = c' * inv(Vnm'*inv(Sxx)*Vnm) * Vnm' * inv(Sxx);
wnm=wnm_conj';
wnmB=conj(B)*wnm; 
% Note: this leads to wnmB*Ynm=conj(y).
% Then plot abs(conj(y))=abs(y).
% This way plot_contour and plot_balloon can be used directly

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


% Additional nulls
th2=70*pi/180; ph2=290*pi/180; % disturbance direction 2
Y2=sh2(N,th2,ph2);
vnm2=B*conj(Y2);

th3=15*pi/180; ph3=310*pi/180; % disturbance direction 2
Y3=sh2(N,th3,ph3);
vnm3=B*conj(Y3);

% constraints
Vnm=[vnm0,vnm1,vnm2,vnm3];
c=[1;0;0;0]

% LCMV with three nulls
wnm_conj = c' * inv(Vnm'*inv(Sxx)*Vnm) * Vnm' * inv(Sxx);
wnm=wnm_conj';
wnmB=conj(B)*wnm; 
% Notes: this computes conj(y) and not y, see above

figure;
plot_contour(wnmB);
hold on; 
plot((180/pi)*ph0,(180/pi)*th0,'w+','MarkerSize',14,'LineWidth',3.0);
plot((180/pi)*ph1,(180/pi)*th1,'k+','MarkerSize',14,'LineWidth',3.0);
plot((180/pi)*ph2,(180/pi)*th2,'k+','MarkerSize',14,'LineWidth',3.0);
plot((180/pi)*ph3,(180/pi)*th3,'k+','MarkerSize',14,'LineWidth',3.0);



figure;
plot_balloon(wnmB,[40,30],0.1); 
axis on; axis tight; axis equal;
xlabel('$x$','FontSize',AxisFontSize,'Interp','Latex');
ylabel('$y$','FontSize',AxisFontSize,'Interp','Latex');
zlabel('$z$','FontSize',AxisFontSize,'Interp','Latex');
set(gca,'FontSize',AxisFontSize);



% Additional distortionless response constraints
th4=55*pi/180; ph4=36*pi/180;
Y4=sh2(N,th4,ph4);
vnm4=B*conj(Y4);
th5=65*pi/180; ph5=36*pi/180;
Y5=sh2(N,th5,ph5);
vnm5=B*conj(Y5);
th6=60*pi/180; ph6=31*pi/180;
Y6=sh2(N,th6,ph6);
vnm6=B*conj(Y6);
th7=60*pi/180; ph7=41*pi/180;
Y7=sh2(N,th7,ph7);
vnm7=B*conj(Y7);

% constraints
Vnm=[vnm0,vnm1,vnm4,vnm5,vnm6,vnm7];
c=[1;0;1;1;1;1]

% LCMV with extanded main lobe 
wnm_conj = c' * inv(Vnm'*inv(Sxx)*Vnm) * Vnm' * inv(Sxx);
wnm=wnm_conj';
wnmB=conj(B)*wnm; 
% Note: this computes conj(y) and not y, see above

figure;
plot_contour(wnmB);
hold on; 
plot((180/pi)*ph0,(180/pi)*th0,'w+','MarkerSize',14,'LineWidth',2.0);
plot((180/pi)*ph4,(180/pi)*th4,'w+','MarkerSize',14,'LineWidth',2.0);
plot((180/pi)*ph5,(180/pi)*th5,'w+','MarkerSize',14,'LineWidth',2.0);
plot((180/pi)*ph6,(180/pi)*th6,'w+','MarkerSize',14,'LineWidth',2.0);
plot((180/pi)*ph7,(180/pi)*th7,'w+','MarkerSize',14,'LineWidth',2.0);
plot((180/pi)*ph1,(180/pi)*th1,'k+','MarkerSize',14,'LineWidth',3.0);

figure;
plot_balloon(wnmB,[40,30],0.1); 
axis on; axis tight; axis equal;
xlabel('$x$','FontSize',AxisFontSize,'Interp','Latex');
ylabel('$y$','FontSize',AxisFontSize,'Interp','Latex');
zlabel('$z$','FontSize',AxisFontSize,'Interp','Latex');
set(gca,'FontSize',AxisFontSize);

% figure(1); print -dpng -r400 ../../figures/chapter07/fig_lcmv1_contour1.png
% figure(2); print -dpng -r600 ../../figures/chapter07/fig_lcmv1_balloon1.png
% figure(3); print -dpng -r400 ../../figures/chapter07/fig_lcmv1_contour2.png
% figure(4); print -dpng -r600 ../../figures/chapter07/fig_lcmv1_balloon2.png
% figure(5); print -dpng -r400 ../../figures/chapter07/fig_lcmv1_contour3.png
% figure(6); print -dpng -r600 ../../figures/chapter07/fig_lcmv1_balloon3.png

