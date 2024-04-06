function [ bobj ] = HeadRotationGenSH( hobj, anm, roomDims,srcAngPos, N, tiledized)
% This function generate BRIR with head rotation for a given HRTF (hobj)
% and plane-waves anm.
%
% Zamir Ben-Hur
% 4.3.2015

if nargin<6
    tiledized=false;
end

%WignerD Matrix
load('Data/WignerDMatrix_diagN=32.mat');
DN=(N+1)^2; % size of the wignerD matrix
D_allAngles = D(:,1:DN);

% SH transform
if strcmp(hobj.dataDomain{2},'SPACE')
    hobj = hobj.toSH(N,'SRC');
else
    warning('hobj is already in the SH domain')
end

% FFT stuff
nFFT = size(anm,2);
if strcmp(hobj.dataDomain{1},'FREQ') && size(hobj.data,2)~=ceil(nFFT/2)+1, hobj=hobj.toTime(); end;
hobj = hobj.toFreq(nFFT);

% Trim negative frequencies
hobj.data = hobj.data(:,1:ceil(nFFT/2)+1,:);
anm = anm(:,1:ceil(nFFT/2)+1,:);
%fVec = fVec(1:(nFFT/2)+1);


% Iterate through head rotation angles
if ~hobj.shutUp, fprintf('Generating BRIRs...\n'); end;
if tiledized
    anm_tilde = anm.';
else
    anm_tilde = anm.'*tildize(N);
end

brirAngles = (pi/180)*(0:1:359);

for jj=1:length(brirAngles)
    
    if ~hobj.shutUp
        if jj==1
            fprintf('   ---> Head rotation = %03d',round(brirAngles(jj)*180/pi)); 
        else
            fprintf('\b\b\b%03d',round(brirAngles(jj)*180/pi));
        end
        
    end
    
    % Rotation matrix
    D = diag(D_allAngles(jj,:));
    
    Hnm_lt_rot=(hobj.data(:,:,1).'*D).';
    Hnm_rt_rot=(hobj.data(:,:,2).'*D).';
    
    % Generate BRIR
    
    plk = sum(anm_tilde.*Hnm_lt_rot.',2).';
    prk = sum(anm_tilde.*Hnm_rt_rot.',2).';
    
    plk(1)=real(plk(1));
    prk(1)=real(prk(1));
    plk(end)=real(plk(end));
    prk(end)=real(prk(end));
    plk=[plk,fliplr(conj(plk(2:end-1)))];
    prk=[prk,fliplr(conj(prk(2:end-1)))];
    arrLt(:,jj)=ifft(plk,'symmetric');
    arrRt(:,jj)=ifft(prk,'symmetric');
end
fprintf('\b\b\bDone!\n')
bobj = BobjGen( hobj, arrLt, arrRt, brirAngles, roomDims,srcAngPos , N);
end

%% Internal functions

function [ bobj ] = BobjGen( hobj, arrLt, arrRt, brirAngles, roomDims,srcAngPos , N)
% Generate the earo object for the BRIR
%

newSrir = earo();
newSrir.name = sprintf('Modeled BRIR in room of dimensions %s',num2str(roomDims));
newSrir.context = sprintf('Rendered BRIR from ISM Model in SH domain + HRTF Set; N=%d.',N);
newSrir.location = 'Virtual room';
newSrir.date = date;
newSrir.engineer = 'Zamir Ben-Hur';
newSrir.contact = 'mail zami@post.bgu.ac.il';
newSrir.comments = sprintf('HRTF data is based on %s',hobj.name);
newSrir.earoVersion = hobj.earoVersion;
newSrir.type = 'BRIR';
newSrir.fs = double(hobj.fs);
newSrir.nData = 360;
newSrir.capturingSystem = 'binauralt_sim.m';
newSrir.sourceGrid.r = srcAngPos(1);
newSrir.positionReference = 'Head Rotation';
newSrir.micGrid.quadType = hobj.sourceGrid.quadType;
newSrir.scatterer = 1;
newSrir.micGrid.r = 0.0875;
newSrir.micGrid.azimuth = brirAngles;
newSrir.micGrid.quadWeight = ones(1,360)*1/360;
newSrir.micGrid.elevation = ones(1,360)*pi/2; % 90 elevation
newSrir.data(:,:,1) = arrLt.';
newSrir.data(:,:,2) = arrRt.';
newSrir.orderN = N;
newSrir.dataDesc = 'BRIR Data head rotation x nData x receivers (1=left ear)';
newSrir.angles = 'RAD';
newSrir.dataDomain{1} = 'TIME';
newSrir.dataDomain{2} = 'SPACE';

bobj=newSrir;
end

function [ Perm ] = tildize( N )
%A_TILD Summary of this function goes here
%   Detailed explanation goes here
Perm=(-1).^(2:(N+1)^2+1);
Perm=diag(Perm);
for n=0:N
    Perm(n^2+1:n^2+2*n+1,n^2+1:n^2+2*n+1)=fliplr(Perm(n^2+1:n^2+2*n+1,n^2+1:n^2+2*n+1));
end
end
