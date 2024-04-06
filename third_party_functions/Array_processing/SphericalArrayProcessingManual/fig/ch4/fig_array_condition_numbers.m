%FIG_ARRAY_CONDITION_NUMBERS generates figures 4.7-4.11, 
% comparing the condition number of matrix B 
% for various array designs, and illustrating 
% the position of samples.
%
% Remark: uses data in fig_array_condition_numbers.mat
% 
% Fundmentals of Spherical Array Processing
% Boaz Rafaely, 2017.

close all;
clear all;

path(path,'../../math');
path(path,'../../plot');

AxisFontSize=14;

clear path;

Na=3; % array order
N=6; % order of sampling scheme
alpha=1.3; % radii ratio
r2=1; % normalized larger radius
r1=r2/alpha; % smaller radius

Nk=500; % wavenumber points
K=linspace(0.1,N,Nk); % wavenumber vector
K=[K,pi+1e-4,4.49341+1e-4,5.7635+1e-4]; % add zeros of b0, b1, b2
K=sort(K);
Nk=length(K);
kr=K*r2;

open_sphere=0; % open sphere
rigid_sphere=1; % rigid sphere
card_mic=2;     % open sphere with cardiod microphones

% Pre-designed sampling points and weights
[a_e,th_e,ph_e]=equiangle_sampling(N);
[a_u,th_u,ph_u]=uniform_sampling(N);
[a_g,th_g,ph_g]=gaussian_sampling(N);
Y_e=sh2(Na,th_e,ph_e);
Y_u=sh2(Na,th_u,ph_u);
Y_g=sh2(Na,th_g,ph_g);
Q=length(th_u);

% single sphere with an additional point at center
r_x=[r2*ones(1,Q),0];
th_ux=[th_u,0];
ph_ux=[ph_u,0];
Y_ux=sh2(Na,th_ux,ph_ux);

% Genetic Algorithm search experiment
% Unmark to run an example of the GA optimization
% V1=ones(1,Q);
% LB=[0*V1];
% UB=[r2*V1];
% options=gaoptimset('InitialPopulation',r2*rand(size(V1)),'Generations',100,'StallTimeLimit',Inf,'TolFun',1e-15,'Display','iter');
% rga = ga(@JkB,Q,[],[],[],[],LB,UB,[],options);
load fig_array_condition_numbers; rga=rga3new; % Load an already optimized data


% Create matrix A (B in paper) at each frequency
for nk=1:Nk;

    k=K(nk);
        
    % single rigid sphere with equiangle sampling
    B=BnMat(Na,k*r2,k*r2,rigid_sphere);
    A_er=diag(B)*Y_e;
    
    % single rigid sphere with gaussian sampling
    A_gr=diag(B)*Y_g;
    
    % single rigid sphere with uniform sampling
    A_ur=diag(B)*Y_u;
    
    % single open sphere with uniform sampling
    B=BnMat(Na,k*r2,k*r2,open_sphere);
    A_uo=diag(B)*Y_u;
    
    %if nk==104, figure, plot_matrix(abs(A_uo)); end;
    if nk==259, A_0mic=A_uo; end;
    
    % single open sphere with uniform sampling, additional point at center
    B=BnMat(Na,k*r_x,k*r_x,open_sphere);
    A_ux=B.'.*Y_ux;
    
    
   % single open sphere, cardioid mic, uniform sampling
   B=BnMat(Na,k*r2,k*r2,card_mic);
   A_uc=diag(B)*Y_u;

   % dual open sphere, uniform sampling, max selection
   Br1=BnMat(Na,k*r1,k*r1,open_sphere);
   Br2=BnMat(Na,k*r2,k*r2,open_sphere);
   B=max(Br1,Br2);
   A_du=diag(B)*Y_u;
    
   % dual open sphere, uniform sampling, no selection
   A_dn=[A_du diag(Br2)*Y_u];
   
  
   % open shell with linearly varying radius r1 to r2
   r=linspace(r1,r2,Q);
   B=BnMat(Na,k*r,k*r,open_sphere);
   A_s1=B.'.*Y_u;
   

   % GA-optimized radii over a ball 0 to r2
   B=BnMat(Na,k*rga,k*rga,open_sphere);
   A_ga=B.'.*Y_u;
      
   
    % condition numbers
    k_er(nk)=cond(A_er);
    k_gr(nk)=cond(A_gr);
    k_ur(nk)=cond(A_ur);
    
    k_uo(nk)=cond(A_uo);
    k_ux(nk)=cond(A_ux);
    
    k_uc(nk)=cond(A_uc);
    k_du(nk)=cond(A_du);
    k_dn(nk)=cond(A_dn);
    
    k_s1(nk)=cond(A_s1);
    k_ga(nk)=cond(A_ga);       

end;
    

figure;
semilogy(kr,k_er,'-','LineWidth',1,'Color',[0 0 0]); hold on;
semilogy(kr,k_gr,'--','LineWidth',1.5,'Color',[0 0 0.5]);
semilogy(kr,k_ur,'-','LineWidth',2.5,'Color',[0 0 0.5]);
hold off;
axis([kr(1) kr(end) 1 10000]);
legend('equiangle','Gaussian','uniform');
xlabel('$kr$','FontSize',AxisFontSize,'Interp','Latex');
ylabel('Condition number','FontSize',AxisFontSize,'Interp','Latex');
set(gca,'FontSize',AxisFontSize);

figure;
semilogy(kr,k_ux,'--','LineWidth',1.5,'Color',[0 0 0.5]); hold on;
semilogy(kr,k_uo,'-','LineWidth',1,'Color',[0 0 0]);
semilogy(kr,k_ur,'-','LineWidth',2.5,'Color',[0 0 0.5]);
hold off;
axis([kr(1) kr(end) 1 10000]);
legend('open + origin','open','rigid','Location','NorthWest'); % corrected
xlabel('$kr$','FontSize',AxisFontSize,'Interp','Latex');
ylabel('Condition number','FontSize',AxisFontSize,'Interp','Latex');
set(gca,'FontSize',AxisFontSize);

figure;
semilogy(kr,k_ur,'-','LineWidth',1,'Color',[0 0 0]); hold on;
semilogy(kr,k_uc,'--','LineWidth',1,'Color',[0 0 0]);
semilogy(kr,k_du,'--','LineWidth',2,'Color',[0 0 0.5]);
semilogy(kr,k_dn,'-','LineWidth',2,'Color',[0 0 0.5]);
hold off;
axis([kr(1) kr(end) 1 10000]);
legend('rigid','cardioid','dual-max','dual-both','Location','NorthEast');
xlabel('$kr$','FontSize',AxisFontSize,'Interp','Latex');
ylabel('Condition number','FontSize',AxisFontSize,'Interp','Latex');
set(gca,'FontSize',AxisFontSize);

figure;
semilogy(kr,k_ur,'-','LineWidth',1,'Color',[0 0 0]); hold on;
semilogy(kr,k_s1,'-','LineWidth',2,'Color',[0 0 0.5]);
semilogy(kr,k_ga,'--','LineWidth',1.5,'Color',[0 0 0.5]);
hold off;
axis([kr(1) kr(end) 1 10000]);
legend('rigid','shell - uni','shell - opt','Location','NorthEast');
xlabel('$kr$','FontSize',AxisFontSize,'Interp','Latex');
ylabel('Condition number','FontSize',AxisFontSize,'Interp','Latex');
set(gca,'FontSize',AxisFontSize);

figure;
h=polarplot(th_u,rga,'bx');
set(h,'LineWidth',2,'Color','blue','MarkerSize',10);
set(gca,'FontSize',AxisFontSize);

figure;
h=polarplot(ph_u,rga,'bx');
set(h,'LineWidth',2,'Color','blue','MarkerSize',10);
set(gca,'FontSize',AxisFontSize);

% figure(1); print -dpng -r400 ../../figures/chapter04/fig_k_rigid.png;
% figure(2); print -dpng -r400 ../../figures/chapter04/fig_k_open.png;
% figure(3); print -dpng -r400 ../../figures/chapter04/fig_k_card_dual.png;
% figure(4); print -dpng -r400 ../../figures/chapter04/fig_k_shell.png;
% figure(5); print -dpng -r400 ../../figures/chapter04/fig_k_rGA_th.png;
% figure(6); print -dpng -r400 ../../figures/chapter04/fig_k_rGA_ph.png;




function J=JkB(r)
%JkBOobjective function for optimization.
% Returns the maximum condition number.

Na=3; % array order
N=6; % order of sampling scheme
alpha=1.3; % radii ratio
r2=1; % normalized larger radius
r1=r2/alpha; % smaller radius
sphere=0; % open sphere
c=343; % speed of sound in m/s

Nk=50; % number of wavenumbers, can be changed 
kmax=N/r2;
kmin=Na/r2;
krange=linspace(kmin,kmax,Nk);

[a,th,ph]=uniform_sampling(N);
Q=length(th);

for nk=1:length(krange),
    k=krange(nk);
    Y=sh2(Na,th,ph);
    Br=BnMat(Na,k*r,k*r,sphere);
    BY=Br.'.*Y;
    con(nk)=cond(BY);
    
end;

J=max(con);

end