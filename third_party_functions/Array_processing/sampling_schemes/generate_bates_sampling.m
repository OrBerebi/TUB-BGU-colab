function [th,ph,L] = generate_bates_sampling(num_points, minTheta, maxTheta, M, plot)
%
% Bates-2015 "Novel Sampling Scheme on the Sphere for Head-Related Transfer FunctionMeasurements"
%
% Inputs:
%         num_points :          wanted number of points (the grid will be as close as possible to this number (the total number of points need to be square rooted)
%         minTheta   :          minimum elevation for the HRTF measurement (in degree, between 0 to 180)
%         maxTheta   :          maximum elevation for the HRTF measurement (in degree, between 0 to 180)
%         M          :          resolution for elevation angles searching  (in degree, for example M=0.1)
%         plot       :          flag for plotting the sampling point
%
% Outputs:
%         th         :          vector of elevation angles (in radians, 0 to pi)
%         ph         :          vector of azimuth angles   (in radians, 0 to 2pi)
%         L          :          maximun spherical harmonics order for the scheme
% Written by Zamir Ben-Hur 29.1.2017
% Edited By Itai Ifergan 23.12.2018 (Changed SH functions and elevation angles to meet ACL convention)
%%%

if nargin<5
    plot = 0; % default dont plot
end

L = ceil(sqrt(num_points));

theta_m = ((minTheta:M:maxTheta)*pi/180).';

% theta_m = flipud(theta_m);
theta = zeros(L,1);


theta(1) = theta_m(1);
theta_m = theta_m(2:end); % taking out zero
% chosing theta
for z = 1 : L-1
    % chosing theta from theta_m which gives the minimum condition number of P
    cond_p = zeros(1,length(theta_m));
    for j = 1:length(theta_m)
        P = [];
        p = shin(L,[theta(1:z); theta_m(j)],zeros(size(theta(1:z),1)+1,1));
        for l = 0:1:z
            plm = p(l^2 + l + 1,:).';
            P = [P plm];
        end
        cond_p(j) = cond(P);
    end
    [~,i] = min(cond_p);
    theta(z+1) = theta_m(i);
    % taking the chosen theta out
    if i==1
        theta_m = theta_m(2:end);
    elseif i==length(theta_m)
        theta_m = theta_m(1:end-1);
    else
        theta_m = [theta_m(1:i-1); theta_m(i+1:end) ];
    end
end

theta_new = zeros(L,1);
[i,~,theta_new(end)] = findnearest(pi/2,theta); % find the closest point to the equator
if i==1
    theta = theta(2:end);
elseif i==length(theta)
    theta = theta(1:end-1);
else
    theta = [theta(1:i-1); theta(i+1:end) ];
end
% reorder vector theta by condition number minimization method (Khalid 2014 �An optimaldimensionalitysampling scheme on the sphere with fast sphericalharmonic transforms�)
for m = L-2:-1:0
    cond_p = zeros(1,length(theta));
    for j = 1:length(theta)
        P = [];
        p = shin(L,[theta(j); theta_new(m+2:end)],zeros(size(theta_new(m+2:end),1)+1,1));
        for l = m:1:L
            plm = p(l^2 + l + m + 1,:).';
            P = [P plm];
        end
        cond_p(j) = cond(P);
    end
    [~,i] = min(cond_p);
    theta_new(m+1) = theta(i);
    % taking the chosen theta out
    if i==1
        theta = theta(2:end);
    elseif i==length(theta)
        theta = theta(1:end-1);
    else
        theta = [theta(1:i-1); theta(i+1:end) ];
    end
end

theta = theta_new;
phi = cell(size(theta));
num_points = 0;
for n = 0:length(theta)-1
    delta_n = 2*pi/(2*n+1);
    phi{n+1} = 0:delta_n:2*n*delta_n;
    num_points = num_points + length(phi{n+1});
end

th = zeros(num_points,1);
ph = zeros(num_points,1);
k=1;
for i = 1:length(theta)
    for j = 1:length(phi{i})
        th(k) = theta(i);
        ph(k) = phi{i}(j);
        k=k+1;
    end
end
L = L-1;
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
end

function Y = shin(N,theta,phi)
% function Y = sh2(N,theta,phi);
%
% N maximum order for n
% theta are the angles theta for the entire location vector
% phi are the angle phi for the entire location vector
% Y is (N+1)^2 by length(theta)

% complex constant
j=sqrt(-1);

% ensure sizes are correct
if length(theta) ~= length(phi)
    fprintf('Lengths of theta and phi must be equal!');
    return;
end
L=numel(theta);
theta=reshape(theta,1,L); % make a row vector
phi=reshape(phi,1,L);       % make a row vector


%n=0
Y=sqrt(1/(4*pi))*ones(1,L);

for n=1:N
    
    % positive m
    a=zeros(1,n);
    for m=0:n
        a(m+1) = sqrt( ((2*n+1)/(4*pi)) * factorial(n-m) / factorial(n+m) );
    end
    a=reshape(a,n+1,1);
    m=(0:n)';
    Y1 = (a*ones(1,L)) .* legendre(n,cos(theta)) .* exp(j*m*phi) ;
    
    % negative m
    m=(-n:-1)';
    Y2 = (((-1).^m)*ones(1,L)) .* conj(Y1(end:-1:2,:));
    
    % append to Y
    Y = [Y; Y2; Y1];
    
end
end


