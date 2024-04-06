function [anm, fVec] = simRoomAnm( N, roomDims, recPos,Orientation,srcPos, srcType, freqAbs, alphaAbs, fs,Anechoic,timeDomainOutput)
    if nargin<11
        timeDomainOutput = false;
    end
%     % collect variables
%     srcAz = srcAngPos(3);
%     srcEl = srcAngPos(2);
%     srcDist = srcAngPos(1);

    % Set-up model
    absMatrix = repmat(alphaAbs, 6, 1);
    if isempty(freqAbs) || isempty(absMatrix)
        Room=SetupRoom('Dim',roomDims);
    else
        Room=SetupRoom('Dim',roomDims,'Freq',freqAbs,'Absorption',absMatrix);
    end
    if (size(recPos,1)==1)
%         Receivers = AddReceiver(Receivers,  'Type',         'gain',     ...
%                                     	'Location',     [2,2,2],    ...
%                                      	'Orientation',  [0,0,0],    ...
%                                    	'Response',     ones(300,3),...
%                                       'Direction',    [az el]     ...
%                            );
        %Receivers=AddReceiver('Location',recPos,'Orientation', [Orientation,0,0],'Type','sphharm','MaxOrder',N,'NFComp',true);
        Receivers=AddReceiver('Location',recPos,'Type','sphharm','MaxOrder',N,'NFComp',true);
    else
        %Receivers=AddReceiver([],'Location',recPos(1,:),'Orientation', [Orientation,0,0],'Type','sphharm','MaxOrder',N,'NFComp',true);
        %Receivers=AddReceiver(Receivers,'Location',recPos(2,:),'Orientation', [Orientation,0,0],'Type','sphharm','MaxOrder',N,'NFComp',true);
        Receivers=AddReceiver([],'Location',recPos(1,:),'Type','sphharm','MaxOrder',N,'NFComp',true);
        Receivers=AddReceiver(Receivers,'Location',recPos(2,:),'Type','sphharm','MaxOrder',N,'NFComp',true);
    end
    %[SrcX SrcY SrcZ] = sph2cart(srcAz,pi/2-srcEl,srcDist);
    
    %Sources=AddSource([],'Location',[SrcX, SrcY, SrcZ] + recPos,'Type',srcType,'Orientation',[180*(pi+srcAz)/pi,0,0]);
    Sources=AddSource([],'Location',srcPos,'Type',srcType);
     % Execute and convert SH types    
  %  simord = [8 8 8];
    if ~Anechoic, 
        Options=MCRoomSimOptions('SimDirect',true,'SimSpec',true,'SimDiff',true,'Fs',double(fs));
    else
        Options=MCRoomSimOptions('SimDirect',true,'SimSpec',false,'SimDiff',false,'Fs',double(fs));
    end
    %if ~hobj.shutUp, PlotSimSetup(Sources,Receivers,Room); end;

    % Execute and convert SH types
    irArray=RunMCRoomSim(Sources,Receivers,Room,Options);
    P=MCRoomPerm(N)';
    C=SHc2r(N)';
    if iscell(irArray)
        anm_L=C*P*irArray{1}.';
        anm_R=C*P*irArray{2}.';
        anm = [anm_L;anm_R];
        
    else
        anm=C*P*irArray.';
    end
    
    
    if ~timeDomainOutput
        % FFT stuff
        nFFT = 2^nextpow2(size(anm,2));
        anm = fft(anm, nFFT, 2);
       % nFFT = size(anm,2);
        fVec = linspace(0,double(fs),nFFT);
    else
        fVec = [];
    end
end

%% Internal functions
function perm=miniSHc2r(n)

    % a help function for SHc2r, permuting for each given n.

    perm = zeros((2*n+1));
    sizeP = size(perm,1);
    perm((floor(sizeP/2)+1),(floor(sizeP/2)+1)) = 1;
    for ii= 1:(floor(sizeP/2))
        perm((floor(sizeP/2)+1+ii),(floor(sizeP/2)+1+ii)) = 1/sqrt(2)*(-1)^ii;%*(-1)^ii;
        perm((floor(sizeP/2)+1+ii),(floor(sizeP/2)+1-ii)) = 1/sqrt(2);
        perm((floor(sizeP/2)+1-ii),(floor(sizeP/2)+1-ii)) = -1/(sqrt(2)*1j);%*(-1)^ii;
        perm((floor(sizeP/2)+1-ii),(floor(sizeP/2)+1+ii)) = +1/(sqrt(2)*1j)*(-1)^ii;
    end
end

function Perm=SHc2r(Nmax)

    % this code forms a permute matrix from the Normalized Complex Spherical Harmonics to
    % the Normalized Real Spherical Harmonics
    % Perm matrix hold the relation- Ynm_{Real} = Perm x Ynm_{Complex}

    Perm = zeros((Nmax+1)^2);
    sizeP = size(Perm,1);
    ind = 0;
    for n= 0:Nmax

        Perm((ind+1):(ind+1+(2*n+1)-1),(ind+1):(ind+1+(2*n+1)-1)) = miniSHc2r(n);
        ind = ind + (2*n +1);
    end

    Perm=conj(Perm);
   
end

function perm=miniMCRoomPerm(n)

    % a help function for MCRooPerm, permuting for each given n.

    perm = zeros((2*n+1));
    sizeP = size(perm,1);
    perm((floor(sizeP/2)+1),(2*n+1)) = 1;
    for ii= 1:(floor(sizeP/2))
        perm((floor(sizeP/2)+1-ii),(2*n+1) - 2*ii +1 ) = 1;
        perm((floor(sizeP/2)+1+ii),(2*n+1) - 2*ii ) = 1;
    end
end

function Perm=MCRoomPerm(Nmax)

    % this code forms a permute matrix that orders the coefficients that we use
    % to the order MCRoomSim does, following GenSHIndices.m The following does
    % so by C_{MCRoomSIM convention} = Perm x C_{our convention}

    Perm = zeros((Nmax+1)^2);
    sizeP = size(Perm,1);
    ind = 0;
    for n= 0:Nmax

        Perm((ind+1):(ind+1  +(2*n+1) - 1     ),(ind+1):(ind+1  +(2*n+1) - 1     )) = miniMCRoomPerm(n);
        ind = ind + (2*n +1);
    end

    Perm = inv(Perm); 
end

function [ Perm ] = tildize( N )
%A_TILD Summary of this function goes here
%   Detailed explanation goes here
    Perm=(-1).^(2:(N+1)^2+1);
    Perm=diag(Perm);
    for n=0:N;
        Perm(n^2+1:n^2+2*n+1,n^2+1:n^2+2*n+1)=fliplr(Perm(n^2+1:n^2+2*n+1,n^2+1:n^2+2*n+1));
    end
end