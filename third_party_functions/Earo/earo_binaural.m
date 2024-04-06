function out = binaural_sh(eobj, hobj, N, a_max)
% function out = binaural_sh(eobj, hobj, N, [a_max])
%
% Render a single binaural impulse response from EARO data
%
%   eobj        - Input EARO object for the room array data
%   hobj        - Input EARO object for the HRTF data
%   N           - Processing order
%   [a_max]     - (optional) radial filter soft-limiting factor in dB, 
%                               {default = inf = no limiting}
%
%   Output:
%   out         - length(f) x 2 array of audio data representing the
%                 rendered BRIR, (:,1) is left ear.
%
% August 2014, Jonathan Sheaffer
% Part of the EARO Beamforming toolbox
%

    fftPad = 1.5;

    if nargin<4, a_max=inf; end;
    
    % Transform to FREQ, if needed
    if strcmp(eobj.dataDomain{1},'TIME')
        if ~eobj.shutUp, fprintf('Transforming to frequency domain...\n'); end;
        nFFT = size(eobj.data,2)*fftPad;
        if mod(nFFT,2), nFFT=nFFT+1; end;   % make sure it is even
        [eobj,fVec] = eobj.toFreq(nFFT);
        fVec=double(fVec);
        if strcmp(hobj.dataDomain{1},'FREQ'), hobj=hobj.toTime(); end;
        hobj = hobj.toFreq(nFFT);   

        % Trim negative frequencies

           eobj.data = eobj.data(:,1:(nFFT/2)+1,:);
           hobj.data = hobj.data(:,1:(nFFT/2)+1,:); 
           fVec = fVec(1:(nFFT/2)+1);
              
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
    anm = earo_pwd(eobj,N,a_max);

    % Compute anm_tilde
    anm_tilde = PWDconj(anm);
    % Y=shMatrix(N,eobj.micGrid.elevation,eobj.micGrid.azimuth);
    % anm_tilde=conj(anm*Y)*pinv(Y);

    % Generate BRIR
    if ~eobj.shutUp, fprintf('Generating BRIR...\n'); end;
    for ii=1:size(anm_tilde,1)
        plk(ii)=anm_tilde(ii,:)*hobj.data(:,ii,1);
        prk(ii)=anm_tilde(ii,:)*hobj.data(:,ii,2);
    end
    plk(1)=real(plk(1));
    prk(1)=real(prk(1));
    plk(end)=real(plk(end));
    prk(end)=real(prk(end));
    plk=[plk,fliplr(conj(plk(2:end-1)))];
    prk=[prk,fliplr(conj(prk(2:end-1)))];
    brir_lt=ifft(plk,'symmetric'); 
    brir_rt=ifft(prk,'symmetric'); 

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
