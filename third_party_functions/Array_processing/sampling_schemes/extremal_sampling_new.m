%EXTREMAL_SAMPLING This function loads an extremal sampling scheme that is
%sufficient for SH representation order N. what is speccial about this
%sampling scheme is that it follows exactly the relation (N+1)^2=Q
%% INPUT
% N             : maximum SH that can accurately represented with this
% type          : type of points set- 'me': minimum energy. 'md': maximum determiant (see details in https://web.maths.unsw.edu.au/~rsw/Sphere/)
% sampling scheme
%% OUTPUT
% r,th,ph        : radius/ elevation and azimuth of each sampling point
%w              : quadrature weights of each sampling point
%% NOTES
% based on Sloan, Ian H., and Robert S. Womersley. "Extremal systems of points and numerical integration on the sphere." Advances in Computational Mathematics 21.1-2 (2004): 107-125.
%
% Edited by Zamir Ben-Hur 18.2.19 - 

function [r,th,ph,w]=extremal_sampling_new(N,type)
if nargin<1
    warning('No order selected. using order N=15')
    N=15;
end
if nargin<2
    type = 'me';
end
if (N<1 || N>40)
    error(['Order N=' num2str(N) ' is not available. N must be betwwen 1 and 40.'])
end

[dir_path,~,~] = fileparts(which('extremal_sampling'));
dir_path=[dir_path,'/Extremal/Extremal system - New2007/' type '/'];
switch type
    case 'me'
        fid = fopen([dir_path, type num2str(N,'%02d') '.' num2str((N+1)^2,'%04d')], 'r');
    case 'md'
       fid = fopen([dir_path, type num2str(N,'%03d') '.' num2str((N+1)^2,'%05d')], 'r');
    otherwise
        error(['No such type ''' type '''!'])
end
if fid==-1
    error(['Can''t read file: ' dir_path,'me' num2str(N,'%03d') '.' num2str((N+1)^2,'%05d')])
end

a = fscanf(fid, '%g %g %g %g', [4 inf]);    % It has two rows now.
a = a';

%transforms into earo soherical coordinates
[ph,th,r] = cart2sph(a(:,1),a(:,2),a(:,3));
ph(ph<0)=ph(ph<0)+2*pi;
th=pi/2-th;
w=a(:,4);

%makes row vector
ph=ph(:).';
th=th(:).';
r=r(:).';
w=w(:).';
fclose(fid);
