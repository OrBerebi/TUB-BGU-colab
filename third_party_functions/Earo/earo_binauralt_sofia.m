function bobj = earo_binauralt_sofia(eobj, hobj, N, a_max)
% function out = binaural_sh(eobj, hobj, N, [a_max])
%
% Render a catalog of binaural impulse responses for head-tracked binaural
%
%   eobj        - Input EARO object for the room array data
%   hobj        - Input EARO object for the HRTF data
%   N           - Processing order
%   [a_max]     - (optional) radial filter soft-limiting factor in dB, 
%                               {default = inf = no limiting}
%
%   Output:
%   bobj        - Output EARO object for BRIRs with head rotations
%
% August 2014, Jonathan Sheaffer
% Part of the EARO Beamforming toolbox
%

    if nargin<4, a_max=inf; end;
    
    % Transform to FREQ, if needed
    if strcmp(eobj.dataDomain{1},'TIME')
        if ~eobj.shutUp, fprintf('Transforming to frequency domain...\n'); end;
        [eobj,fVec] = eobj.toFreq();
        nFFT = size(fVec,2);
        fVec=double(fVec);
        if strcmp(hobj.dataDomain{1},'FREQ'), hobj=hobj.toTime(); end;
        hobj = hobj.toFreq(nFFT);   

        % Trim negative frequencies

        eobj.data = eobj.data(:,1:(nFFT/2)+1,:);     
        fVec = fVec(1:(nFFT/2)+1);
        hobj.data = hobj.data(:,1:(nFFT/2)+1,:);
        if size(fVec,2)==nFFT, fVec = fVec(1:(nFFT/2)+1); end;
  
    end

    % Transform to SH domain
    if strcmp(eobj.dataDomain{2},'SH') || strcmp(eobj.dataDomain{2},'SH')
        fprintf('Please ensure that both array data and HRTFs are in the space domain.\nSomething is not right.');
    end
    if ~eobj.shutUp, fprintf('Transforming to SH domain, order N = %d:\n',N); end;
    eobj = eobj.toSH(N,'MIC');
    hobj = hobj.toSH(N,'SRC');

    % Do PWD
    if ~eobj.shutUp, fprintf('Performing plane wave decomposition...\n'); end;
    Pnm=squeeze(eobj.data).';
    [compositeGrid, Npwd] = getCompositeGrid(Pnm, true);
    kr=2*pi*fVec*eobj.micGrid.r(1)/343.5;
    dn   = sofia_mf(N , double(kr), 1, a_max);
    A    = sofia_pdc(Npwd, compositeGrid(:,1:2), Pnm, dn);  % 74 x 8501
    Y = shMatrix(N,compositeGrid(:,2),compositeGrid(:,1));
    anm_tilde = conj(A).'*pinv(Y);
    
    % Iterate through head rotation angles
    if ~eobj.shutUp, fprintf('Generating BRIRs...\n'); end;
    brirAngles=(pi/180)*(0:1:359);
    for jj=1:length(brirAngles)

        if ~eobj.shutUp       
            fprintf('   --->Head rotation = %d \r',brirAngles(jj)*180/pi); 
            
        end
        
        % Rotation matrix
        D=WignerDM(N,0,0,brirAngles(jj));
        Hnm_lt_rot=(hobj.data(:,:,1).'*D).';
        Hnm_rt_rot=(hobj.data(:,:,2).'*D).';   
        
        % Generate BRIR
        clear plk prk;
        for ii=1:size(anm_tilde,1)
            plk(ii)=anm_tilde(ii,:)*Hnm_lt_rot(:,ii);
            prk(ii)=anm_tilde(ii,:)*Hnm_rt_rot(:,ii);
        end
        plk(1)=real(plk(1));
        prk(1)=real(prk(1));
        plk(end)=real(plk(end));
        prk(end)=real(prk(end));
        plk=[plk,fliplr(conj(plk(2:end-1)))];
        prk=[prk,fliplr(conj(prk(2:end-1)))];
        arrLt(:,jj)=ifft(plk,'symmetric'); 
        arrRt(:,jj)=ifft(prk,'symmetric');         
    end
    
    newSrir = earo();
    newSrir.name = eobj.name;
    newSrir.context = sprintf('Rendered BRIR from Spherical Array + HRTF Set; N=%d.',N);
    newSrir.location = eobj.location;
    newSrir.date = date;
    newSrir.engineer = 'Jonathan Sheaffer';
    newSrir.contact = 'mail sheaffer@ee.bgu.ac.il';
    newSrir.comments = sprintf('HRTF data is based on %s',hobj.name);
    newSrir.earoVersion = eobj.earoVersion;
    newSrir.type = 'BRIR';
    newSrir.fs = eobj.fs;
    newSrir.nData = 360;
    newSrir.capturingSystem = 'binauralt_sh.m';
    newSrir.sourceGrid.r = eobj.sourceGrid.r;
    newSrir.positionReference = 'Head Rotation';
    newSrir.micGrid.quadType = 'Gauss-Leg. 360SP (1E/360A)';
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

% Internal functions

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

function [f2]=PWDconj(f)

    f=f.';
    N=sqrt(size(f,1))-1;
    for n=0:N
        m=[-n:n]';
        f2(n^2+n+m+1,:)=repmat((-1).^(-m),1,size(f,2)).*(f(n^2+n-m+1,:));
    end

    f2=f2.';

end

function [compositeGrid, Npwd] = getCompositeGrid(Pnm, cGridType)

    Npwd = sqrt(size(Pnm,1))-1;

    if cGridType

        lebNodes = [6, 14, 26, 38, 50, 74, 86, 110, 146, 170, 194, 230, 266, 302, ...
            350, 434, 590, 770, 974, 1202, 1454, 1730, 2030, 2354, 2702, ...
            3074, 3470, 3890, 4334, 4802, 5294, 5810];

        lebOrders  = floor(sqrt(lebNodes/1.3)-1);
        matchOrder = abs(lebOrders-Npwd);
        orderIndex = find(matchOrder==min(matchOrder));

        if ~matchOrder(orderIndex)
            fprintf(['>>> composite grid: Lebedev, ', num2str(lebNodes(orderIndex)),' Nodes\n\n']);
            compositeGrid = sofia_lebedev(lebNodes(orderIndex),0);
            return
        else
            warning('No matching Lebedev composit grid available: Switched to Gauss grid.\n');
        end
    end

    fprintf(['>>> composite grid: Gauss, ', num2str(2*(Npwd+1)^2),' Nodes\n\n']);
    compositeGrid = sofia_gauss(2*(Npwd+1), (Npwd+1),0);
end