function [a,th,ph,max_SHorder]=generate_gaussian_sampling_numOfPoints(num_points, minTheta, maxTheta, plot)

% [a,th,ph]=generate_gaussian_sampling(N);
%
% a - weights
% th - theta angles for all points
% ph - phi angles for all points
% Total no. of points is length(a)=length(th)=length(ph)
%
% Gaussina sampling for order N
%
% Boaz Rafaely 12 October 2006
% updated by Zamir Ben-Hur 11 November 2016 
if nargin<4
    plot = 0; % default dont plot
end
if ~exist('maxTheta','var')
    maxTheta = 90;
end
if ~exist('minTheta','var')
    minTheta = -90;
end

N = ceil(sqrt(num_points/2)-1);
num_found_points = 0;

while (num_found_points < num_points)
    Lg=N+1;
    Leq=2*(N+1);
    
    % Equiangle samples over phi
    phi=0:2*pi/Leq:(Leq-1)*2*pi/Leq;
    
    
    % Gaussian samples and nodes over theta
    
    [x,w]=lgwt(N+1,-1,1);
    theta=acos([x; -x(end-1:-1:1)]);
    aa=(pi/(N+1))*[w; w(end-1:-1:1)];
    
    count=0;
    th = zeros(1,Lg*Leq);
    ph = zeros(1,Lg*Leq);
    a = zeros(1,Lg*Leq);
    for j=1:Lg
        for k=1:Leq
            count=count+1;
            th(count) = theta(j);
            ph(count) = phi(k);
            a(count)  = aa(j);
        end
    end
    
    th = pi/2 - th ;
    ind_th = find(th>=minTheta*pi/180 & th<=maxTheta*pi/180);
    num_found_points = length(ind_th);
    N = N+1;
end
th = th(ind_th).';
ph = ph(ind_th).';
a = a(ind_th).';
max_SHorder = floor(sqrt(length(th)/2)-1);

if plot    
    plotSampleDirections([ph th],[],' ',1,1,8); colorbar off
    colormap Gray;
    hold on;
    sphere;
    axis equal;
    view([40 -10]); 
    rotate3d on;
    light;
    alpha(.8);
    lighting phong;
    hold off;
end
