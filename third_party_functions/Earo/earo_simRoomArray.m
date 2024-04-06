
function [eobj, anm_L,anm_R, fVec] = earo_simRoomArray(N, roomDims, recPos_L,recPos_R, srcAngPos, srcType, freqAbs, alphaAbs, fs, array_radius, sphereType)

    % function [eobj, anm, fVec] = earo_simRigidArray(N, roomDims, recPos, 
    %                                       srcAngPos, srcType, freqAbs, alphaAbs, fs, array_radius)

    % See example2.m for input parameters
    %
    % Returns: 
    %    eobj - EARO object for the simulated array in the room
    %    anm - Plane wave density function directly from the simulation
    %    (i.e. no effects of sampling)
    %    fVec - Frequency vector for anm's
    
    
    fftPad = 1;

    % collect variables
    srcAz = srcAngPos(3);
    srcEl = srcAngPos(2);
    srcDist = srcAngPos(1);

    % Set-up model
    absMatrix = repmat(alphaAbs, 6, 1);
    if isempty(freqAbs) || isempty(absMatrix)
        Room=SetupRoom('Dim',roomDims);
    else
        Room=SetupRoom('Dim',roomDims,'Freq',freqAbs,'Absorption',absMatrix);
    end
    Receivers=AddReceiver('Location',recPos_L,'Type','sphharm','MaxOrder',N,'NFComp',false);
    Receivers=AddReceiver(Receivers,'Location',recPos_R,'Type','sphharm','MaxOrder',N,'NFComp',false);

    %[SrcX ,SrcY ,SrcZ] = sph2cart(srcAz,pi/2-srcEl,srcDist);
    %Sources=AddSource([],'Location',[SrcX, SrcY, SrcZ] + recPos,'Type',srcType,'Orientation',[radtodeg(pi+srcAz),0,0]);
    Sources=AddSource([],'Location',srcAngPos,'Type',srcType,'Orientation',[rad2deg(pi+srcAz),0,0]);

    Options=MCRoomSimOptions('SimDirect',true,'SimSpec',true,'SimDiff',true,'Fs',double(fs));
    PlotSimSetup(Sources,Receivers,Room);

    % Execute and convert SH types
    irArray=RunMCRoomSim(Sources,Receivers,Room,Options);
    P=MCRoomPerm(N)';
    C=SHc2r(N)';
    anm_L=C*P*irArray{1,1}.';
    anm_R=C*P*irArray{2,1}.';
    
    % FFT stuff  
    nFFT = size(anm_L,2)*fftPad;
    fVec = linspace(0,double(fs),nFFT); 
    if mod(nFFT,2), nFFT=nFFT+1; end;   % make sure it is even
    
    anm_L = fft(anm_L, nFFT, 2);
    anm_R = fft(anm_R, nFFT, 2);

    % Trim negative frequencies 
    anm_L = anm_L(:,1:(nFFT/2)+1,:);
    anm_R = anm_R(:,1:(nFFT/2)+1,:);
    fVec = fVec(1:(nFFT/2)+1);
    
    kr = fVec*2*pi*array_radius/343.5;   % Default c for mcroomsim
    kr(1)=kr(2);                         % Remove DC
    
    % Generate Pnm
    bn = radialMatrix(N,kr,sphereType,inf,true,false);
    pnm = anm_L.*bn.';
    
    % Back to space
    [aW,thRec,phRec] = generate_gaussian_sampling(N);
    Y=shMatrix(N,thRec,phRec);
    p = Y.'*pnm;
    
    % Back to time
    p(:,1)=real(p(:,1));
    p(:,end)=real(p(:,end));
    p=[p,fliplr(conj(p(:,2:end-1)))];
    irArray=ifft(p,[],2,'symmetric'); 

    newSrir = earo();
    newSrir.name = sprintf('Modeled array (ISM) in room of dimensions %s',num2str(roomDims));
    newSrir.context = sprintf('Array is Gaussian corresponding to N=%d.',N);
    newSrir.location = 'Virtual room';
    newSrir.date = date;
    newSrir.engineer = 'Someone';
    newSrir.contact = 'mail someone@ee.bgu.ac.il';
    newSrir.earoVersion = 1.0;
    newSrir.type = 'arrayData';
    newSrir.fs = fs;
    newSrir.nData = length(aW);
    newSrir.capturingSystem = 'earo_simRigifArray.m';
    newSrir.sourceGrid.r = srcAngPos(1);

    newSrir.micGrid.quadType = 'Gauss-Legendre';
    newSrir.scatterer = 1;
    newSrir.micGrid.r = array_radius;
    newSrir.micGrid.azimuth = phRec;
    newSrir.micGrid.quadWeight = aW;
    newSrir.micGrid.elevation = thRec; % 90 elevation
  
    newSrir.data(1,:,:) = irArray.';
    newSrir.isRobot = false;

    newSrir.orderN = N;
    newSrir.dataDesc = 'src x nData x receivers';
    newSrir.angles = 'RAD';
    newSrir.dataDomain{1} = 'TIME';
    newSrir.dataDomain{2} = 'SPACE';

    eobj=newSrir;    
    
end

%% Internal functions
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


function DM=WignerDM(N,alpha,beta,gamma)
% build Wigner D matrix (eq. (20), Ben Hagai Nov 2012) 
% first, gamma about z, then beta about y, then gamma about z again (counterclocwise)
DM=zeros((N+1)^2);
    for n=0:N
        for m=-n:n
            DM(n^2+n+m+1,n^2+1:(n+1)^2)=...
                wignerd(n,m,alpha,beta,gamma);
        end
    end
end

function D=wignerd(n,m,alpha,beta,gamma)
% wignerd.m
% ------------
%       Get Wigner-D Coefficients.
%       Based on Rafaely 2008, equations (11-12).
%
% Syntax
% ------------
%     D=wignerd(n,m,m2,alpha,beta,gamma)
%
% Input
% ------------
%     Required
%           n,m - scalar - the spherical indices
%           alpha,beta,gamma - Euler angles
%
% Output
% ------------
%         D - defined as D_{mm'}^n, where m'=-n:n
%
% Created/Modified by
% ------------
%     Ilan Ben Hagai, 1-Nov-2010

%%
    m2=-n:n;
    
    epsilon=1.*(m2>=m) + (-1).^(m2-m).*(m2<m);
    mu=abs(m-m2);
    nu=abs(m+m2);
    s=n-(mu+nu)/2;
    Ps=zeros(1,numel(mu));
    for mIdx=1:numel(mu)
        curMu=mu(mIdx);
        curNu=nu(mIdx);
        curS=s(mIdx);
        polyCoeffs=orth_poly('Jacobi',curS,curMu,curNu);
        Ps(mIdx)= polyval(polyCoeffs,cos(beta));
    end
    
    % calculate the Wigner-d function (eq.12)
    d=epsilon.*sqrt(factorial(s).*factorial(s+mu+nu)./(factorial(s+mu).*factorial(s+nu))).*sin(beta/2).^mu.*cos(beta/2).^nu.*Ps;

    % Calculate the coefficients (eq.11) :
    D=exp(-1i*m*alpha-1i*m2*gamma).*d;
end

function pn=orth_poly(class,n,alpha,beta)
% generates an orthogonal polynomial

    if (nargin<4)||isempty(beta)
        beta=0;
    end
    if (nargin<3)||isempty(alpha);
        alpha=0;
    end

    % initialize (-1)'th and zero'th order polynomials
    pn=[];
    pnp1=1;

    for i=0:n
        pnm1=pn;
        pn=pnp1;
        switch class
            case 'Legendre'
                pnp1=((2*i+1)*[pn,0] - i*[0,0,pnm1])/(i+1);
            case 'Hermite'
                pnp1=2*[pn,0] - 2*i*[0,0,pnm1];
            case 'Laguerre'
                pnp1=((2*i+alpha+1)*[0,pn] -[pn,0] - (i+alpha)*[0,0,pnm1])/(i+1);
            case 'Jacobi'
                if (alpha~=0)||(beta~=0)
                    a1n=2*(i+1)*(i+alpha+beta+1)*(2*i+alpha+beta);
                    a2n=(2*i+alpha+beta+1)*(alpha^2-beta^2);
                    if (2*i+alpha+beta)<=150
                        a3n=gamma(2*i+alpha+beta+3)./gamma(2*i+alpha+beta);
                    else
                        a3n=exp(gammaln(2*i+alpha+beta+3)-gammaln(2*i+alpha+beta));
                    end
                    a4n=2*(i+alpha)*(i+beta)*(2*i+alpha+beta+2);
                    pnp1=(a2n*[0,pn] + a3n*[pn,0] - a4n*[0,0,pnm1])./a1n;
                else
                    pnp1=((2*i+1)*[pn,0] - i*[0,0,pnm1])/(i+1);
                end
        end

    end
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

function [a,th,ph]=generate_gaussian_sampling(N)

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


    Lg=N+1;
    Leq=2*(N+1);

    % Equiangle samples over phi
    phi=0:2*pi/Leq:(Leq-1)*2*pi/Leq;


    % Gaussian samples and nodes over theta

    [x,w]=lgwt(N+1,-1,1);
    theta=acos([x; -x(end-1:-1:1)]);
    aa=(pi/(N+1))*[w; w(end-1:-1:1)];

    count=0;
    for j=1:Lg,
        for k=1:Leq,
            count=count+1;
            th(count) = theta(j);
            ph(count) = phi(k);
            a(count)  = aa(j);
        end;
    end
end

function [x,w]=lgwt(N,a,b)

    % lgwt.m
    %
    % This script is for computing definite integrals using Legendre-Gauss 
    % Quadrature. Computes the Legendre-Gauss nodes and weights  on an interval
    % [a,b] with truncation order N
    %
    % Suppose you have a continuous function f(x) which is defined on [a,b]
    % which you can evaluate at any x in [a,b]. Simply evaluate it at all of
    % the values contained in the x vector to obtain a vector f. Then compute
    % the definite integral using sum(f.*w);
    %
    % Written by Greg von Winckel - 02/25/2004
    N=N-1;
    N1=N+1; N2=N+2;

    xu=linspace(-1,1,N1)';

    % Initial guess
    y=cos((2*(0:N)'+1)*pi/(2*N+2))+(0.27/N1)*sin(pi*xu*N/N2);

    % Legendre-Gauss Vandermonde Matrix
    L=zeros(N1,N2);

    % Derivative of LGVM
    Lp=zeros(N1,N2);

    % Compute the zeros of the N+1 Legendre Polynomial
    % using the recursion relation and the Newton-Raphson method

    y0=2;

    % Iterate until new points are uniformly within epsilon of old points
    while max(abs(y-y0))>eps


        L(:,1)=1;
        Lp(:,1)=0;

        L(:,2)=y;
        Lp(:,2)=1;

        for k=2:N1
            L(:,k+1)=( (2*k-1)*y.*L(:,k)-(k-1)*L(:,k-1) )/k;
        end

        Lp=(N2)*( L(:,N1)-y.*L(:,N2) )./(1-y.^2);   

        y0=y;
        y=y0-L(:,N2)./Lp;

    end

    % Linear map from[-1,1] to [a,b]
    x=(a*(1-y)+b*(1+y))/2;      

    % Compute the weights
    w=(b-a)./((1-y.^2).*Lp.^2)*(N2/N1)^2;
end