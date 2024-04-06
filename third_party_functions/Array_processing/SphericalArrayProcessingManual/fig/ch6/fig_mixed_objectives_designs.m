%FIG_MIXED_OBJECTIVES_DESIGNS generates results for Table 6.2, 
% comparing several array designs with mixed objectives.
%
% Fundmentals of Spherical Array Processing
% Boaz Rafaely, 2017.

close all;
clear all;
clc;

path(path,'../../math');
path(path,'../../plot');

% Parameters: sphere, N, kr, Sigma^2a, Sigma^2s
par=[0,2,2,1,0;
     0,2,2,0,1;
     1,3,3,1,0;
     1,3,3,0,1;
     1,3,3,1,1;
     1,4,2,1,0;
     1,4,2,0,1;
     1,4,2,0.4,1];
    
for i=1:size(par,1),

    sphere=par(i,1);
    N=par(i,2);
    [a,th,ph]=uniform_sampling(N);
    Q(i)=length(a);
    NN=(2*[0:N]+1).';
    vn=NN/(4*pi);
   
    kr=par(i,3);
    bn=BnMat(N,kr,kr,sphere); bn=bn([1:N+1].^2).';
    sigma2a=par(i,4);
    sigma2s=par(i,5);
    
    R = sigma2a*(1/(4*pi))*diag(NN) + sigma2s*(1/Q(i))*diag(NN./abs(bn).^2);
    dn=inv(R)*vn/(vn'*inv(R)*vn);
    DI(i)=(dn'*NN)^2/(dn'*diag(NN)*dn);
    WNG(i)=(Q(i)/(4*pi)^2)*(dn'*NN)^2/(dn'*diag(NN./abs(bn).^2)*dn);
end;

% Table
fprintf('\nSph \tN  \tQ  \tkr \tsig_a  \tsig_s  \tDF \t\tWNG');
for i=1:size(par,1),
 fprintf('\n%1.0f & \t%1.0f & \t%1.0f & \t%1.0f &  \t%1.1f &  \t%1.1f &  \t%05.2f & \t%05.2f',par(i,1),par(i,2),Q(i), par(i,3),par(i,4),par(i,5),DI(i),WNG(i));
end;
fprintf('\n');
