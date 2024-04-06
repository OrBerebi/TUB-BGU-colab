function [ITD,ILD,plk,prk,fVec]=bnLimitedCues(hobj,N,a_max,el0,az0)


% Generate binaural cues for radial function limited field


    fs=hobj.fs;
    padFFT = 4;   % FFT padding factor

    % Transform to FREQ, if needed
    if strcmp(hobj.dataDomain{1},'FREQ'), hobj=hobj.toTime(); end;
    [hobj,fVec] = hobj.toFreq(hobj.taps*padFFT);   
    nFFT = size(fVec,2);
    fVec=double(fVec);

    % Trim negative frequencies
    hobj.data = hobj.data(:,1:(nFFT/2)+1,:); 
    fVec = fVec(1:(nFFT/2)+1);
    kr=2*pi*hobj.micGrid.r(1)*fVec/343.5;

    % Transform to SH domain
    if ~strcmp(hobj.dataDomain{2},'SH')
        if ~hobj.shutUp, fprintf('Transforming to SH domain, order N = %d:\n',N); end;
        hobj = hobj.toSH(N,'SRC');    
    end

    
    % Generate a single, bn-limited plane wave
    if a_max~=inf
        Y=conj(shMatrix(N,el0,az0));
        bn_lim=radialMatrix(N,kr,1,a_max,true,true);
        bn_reg=radialMatrix(N,kr,1,inf,true,true);

        Y2=repmat(Y,1,length(fVec)).';
        anm = Y2.*(bn_reg./bn_lim);
    else
        Y=conj(shMatrix(N,el0,az0));
        anm=repmat(Y,1,length(fVec)).';
    end
    anm_tilde = anm*tildize(N);
    
    % Generate BRIR
    if ~hobj.shutUp, fprintf('Generating BRIR...\n'); end;
    for ii=1:size(anm_tilde,1)
        plk(ii)=anm_tilde(ii,:)*hobj.data(:,ii,1);
        prk(ii)=anm_tilde(ii,:)*hobj.data(:,ii,2);
    end
    plk(1)=real(plk(1));
    prk(1)=real(prk(1));
    plk(end)=real(plk(end));
    prk(end)=real(prk(end));
    
    % Compute cues
    ILD = 20*log10(abs(plk./prk));
    ITD = -1./(2*pi*fVec).*angle(plk./prk);
end   
