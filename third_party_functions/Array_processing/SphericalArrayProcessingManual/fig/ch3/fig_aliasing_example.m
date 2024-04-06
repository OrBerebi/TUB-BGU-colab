%FIG_ALIASING_EXAMPLE generates Figures 3.7, 3.8, showing balloon
% plots to illustrate aliasing.
%
% Fundmentals of Spherical Array Processing
% Boaz Rafaely, 2017.

close all
clear all

path(path,'../../math');
path(path,'../../plot');

q = @(n,m) (n^2+n+m+1);

N=5;

f1nm=zeros((N+1)^2,1);
f1nm(q(0,0))=sqrt(4*pi);
f2nm=zeros((N+1)^2,1);
f2nm(q(5,-5))=0.5*sqrt(1024*pi/693);
f2nm(q(5,5))=-0.5*sqrt(1024*pi/693);

f1hat_nm=zeros((N+1)^2,1);
f1hat_nm(q(0,0))=sqrt(4*pi);
f2hat_nm=zeros((N+1)^2,1);
f2hat_nm(q(3,-3))=0.5*sqrt(64*pi/35);
f2hat_nm(q(3,3))=-0.5*sqrt(64*pi/35);

AX=2*[-1,1,-1,1,-1,1];

Zoom=2.5;

figure(1);

subplot(2,2,[1,2]);
plot_balloon(f1nm+f2nm);
axis(AX); axis square; 
camzoom(Zoom);

subplot(223);
plot_balloon(f1nm);
axis(AX); axis square; 
camzoom(Zoom);

subplot(224);
plot_balloon(f2nm);
axis(AX); axis square; 
camzoom(Zoom);


figure(2);

subplot(2,2,[1,2]);
plot_balloon(f1hat_nm+f2hat_nm);
axis(AX); axis square; 
camzoom(Zoom);

subplot(223);
plot_balloon(f1hat_nm);
axis(AX); axis square; 
camzoom(Zoom);

subplot(224);
plot_balloon(f2hat_nm);
axis(AX); axis square;
camzoom(Zoom);



% Perform Equal-angle sampling and check aliasing
NN=3;
[a3,th3,ph3]=equiangle_sampling(NN);
Y1=sh2(N,th3,ph3);
Y=sh2(NN,th3,ph3);
alpha=sqrt(4*pi);
gamma=0.5*sqrt(1024*pi/693);
f=alpha*Y1(q(0,0),:)+gamma*Y1(q(5,-5),:)-gamma*Y1(q(5,5),:);
fnm=conj(Y)*diag(a3)*f.';
 
% figure(1); print -dpng -r600 ../../figures/chapter03/fig_aliasing_example_f.png
% figure(2); print -dpng -r600 ../../figures/chapter03/fig_aliasing_example_f_hat.png
