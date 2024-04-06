function [earoOut, outFilt] = earo_eqBRIR (eobj, targetN, eqMethod, norml, hobj, azEl)

  useSOFiA = false;                     % change to true if you wish to generate the radial functions using SOFiA toolbox
                                        % instead of EARO's internal func.
  timeDep = true;                       % Time dependence assumption (true = positive, i.e. e^jwt); not available with SOFiA(!)


  % Transform to frequency domain, if needed
  if strcmp(eobj.dataDomain{1},'TIME')
      if ~eobj.shutUp, disp('Transforming to frequency domain...'); end;
      nFFT = 2^nextpow2(round(size(eobj.data,2)*1.03)); % Pad to avoid circular conv
      [eobj,fVec]=eobj.toFreq(nFFT); 
  else
      disp('ERROR! BRIRs must be in the time domain');
  end
  
  if nargin>4
    if strcmp(hobj.dataDomain{1},'TIME')
        if ~hobj.shutUp, disp('Transforming to frequency domain...'); end;
        [hobj,tmp]=hobj.toFreq(nFFT); 
    end
  end

  % Transform to space domain, if needed
  if strcmp(eobj.dataDomain{2},'SH')
      if ~eobj.shutUp, disp('Transforming to space domain...'); end;
      eobj=eobj.toSpace('MIC');
  end
  
  % Compute speed of sound
  if ~isempty(eobj.avgAirTemp) && isnumeric(eobj.avgAirTemp)
    c = 331.3 * sqrt(1+(eobj.avgAirTemp/273.15));
  else
      c = 343.5; % default value;
  end
  
  % Trim negative frequencies
    eobj.data = eobj.data(:,1:(nFFT/2)+1,:);
    if nargin>4, hobj.data = hobj.data(:,1:(nFFT/2)+1,:); end;
    fVec = fVec(1:(nFFT/2)+1); 
    
  % Construct kr vector
  fVec=double(fVec);
  kr=double(fVec*2*pi*eobj.micGrid.r(1)/c);
  
  if strcmp(eqMethod,'spherical')      % Equalize using a spherical head model
      
    if ~eobj.shutUp, fprintf('Designing filter, type %s...\n',eqMethod); end;
    % Design spherical-head filter
    sResp = meanRespN(eobj.orderN, kr, timeDep, useSOFiA);
    tResp = meanRespN(targetN, kr, timeDep, useSOFiA);
    sFilt = repmat(tResp./sResp,1,360).';
    outFilt = tResp./sResp;
    sFilt(isnan(sFilt)) = 1;
    outFilt(isnan(outFilt)) = 1;
    
    % Filter data
    if ~eobj.shutUp, disp('Convolving with data and performing IFFT...'); end;
    plk = eobj.data(:,:,1).* sFilt;
    prk = eobj.data(:,:,2).* sFilt;
      
  elseif strcmp(eqMethod,'realhead')    % Equalize using hrtf-based model
    
    % Source Response
    tmp=hobj.toSH(eobj.orderN,'SRC');
    hrtf_lt = tmp.data(:,:,1);
    hrtf_rt = tmp.data(:,:,2);  
    srcLt = sqrt((1/(4*pi))*sum(abs(hrtf_lt).^2,1)).';
    srcRt = sqrt((1/(4*pi))*sum(abs(hrtf_rt).^2,1)).';
    
    % Target Response
    tmp=hobj.toSH(targetN,'SRC');
    hrtf_lt = tmp.data(:,:,1);
    hrtf_rt = tmp.data(:,:,2);  
    tarLt = sqrt((1/(4*pi))*sum(abs(hrtf_lt).^2,1)).';
    tarRt = sqrt((1/(4*pi))*sum(abs(hrtf_rt).^2,1)).';
    
    % Design filters
    sFiltLt = repmat(tarLt./srcLt,1,360).';
    sFiltRt = repmat(tarRt./srcRt,1,360).';
    
    outFilt(:,1)=tarLt./srcLt;
    outFilt(:,2)=tarRt./srcRt;
    
    % Filter data
    if ~eobj.shutUp, disp('Convolving with data and performing IFFT...'); end;
    plk = eobj.data(:,:,1).* sFiltLt;
    prk = eobj.data(:,:,2).* sFiltRt;  
    
  elseif strcmp(eqMethod,'single')         % Equalize using single HRTF
      
    % Truncate HRTFs
    tmp = hobj.toSH(eobj.orderN, 'SRC');
    hobjLow = tmp.toSpace('SRC');

    tmp = hobj.toSH(targetN,'SRC');
    hobjHi = tmp.toSpace('SRC');

    % Find the azimuth/elevation
    dirIdxHi = hobjHi.closestIrSrc(azEl(2), azEl(1));
    dirIdxLo = hobjLow.closestIrSrc(azEl(2), azEl(1));

    % Design the filter
    outFilt(:,1) = abs(hobjHi.data(dirIdxHi,:,1)) ./ abs(hobjLow.data(dirIdxLo,:,1));
    outFilt(:,2) = abs(hobjHi.data(dirIdxHi,:,2)) ./ abs(hobjLow.data(dirIdxLo,:,2));
    sFiltLt = repmat(outFilt(:,1),1,360).';
    sFiltRt = repmat(outFilt(:,2),1,360).';      
      
    % Filter data
    if ~eobj.shutUp, disp('Convolving with data and performing IFFT...'); end;
    plk = eobj.data(:,:,1).* sFiltLt;
    prk = eobj.data(:,:,2).* sFiltRt;          
      
  end
  
  % Make sure that DC and Nyquist are real
  plk(:,1)=real(plk(:,1));
  plk(:,end)=real(plk(:,end));
  prk(:,1)=real(prk(:,1));
  prk(:,end)=real(prk(:,end));
  
  % Generate negative spectrum and transform to time
  plk = [plk, fliplr(conj(plk(:,2:end-1)))];
  prk = [prk, fliplr(conj(prk(:,2:end-1)))];
  
  pl = ifft(plk,[],2,'symmetric');
  pr = ifft(prk,[],2,'symmetric');
  
  % Normalize
  if norml
      nval=max(max(abs(pl))); 
      pl = pl./nval;
      pr = pr./nval;
  end
  
  % Generate output object
  earoOut = eobj;
  earoOut.context = strcat(eobj.context,sprintf(' Eq to N=%d using method: %s',targetN,eqMethod));
  earoOut.data = [];
  earoOut.data(:,:,1) = pl;
  earoOut.data(:,:,2) = pr;
  earoOut.dataDomain{1} = 'TIME';
  earoOut.dataDomain{2} = 'SPACE';

end

%% Internal functions

function out=meanRespN(N, kr, timeDep, useSOFiA)
% function out=meanRespN(N, kr)
%
% Calculate averaged magnitude pressure at left and right ears
% due to diffraction around a rigid sphere, for order N


  bn=radialMatrix (N,kr,1,inf,timeDep,useSOFiA);
  Yl = shMatrix(N,pi/2,pi/2);
  
  
  tmp=(abs(bn).^2)*(abs(Yl).^2);
  out=sqrt(tmp/(4*pi));
  

end


function [r,c,V] = findnearest(srchvalue,srcharray,bias)

% Usage:
% Find the nearest numerical value in an array to a search value
% All occurances are returned as array subscripts
%
% Output:
%
% For 2D matrix subscripts (r,c) use:
%
%       [r,c] = findnearest(srchvalue,srcharray,gt_or_lt)
%
%
% To also output the found value (V) use:
%
%       [r,c,V] = findnearest(srchvalue,srcharray,gt_or_lt)
%
%
% For single subscript (i) use:
%
%         i   = findnearest(srchvalue,srcharray,gt_or_lt)
% 
%
% Inputs:
%
%    srchvalue = a numerical search value
%    srcharray = the array to be searched
%    bias      = 0 (default) for no bias
%                -1 to bias the output to lower values
%                 1 to bias the search to higher values
%                (in the latter cases if no values are found
%                 an empty array is ouput)
%
%
% By Tom Benson (2002)
% University College London
% t.benson@ucl.ac.uk

    if nargin<2
        error('Need two inputs: Search value and search array')
    elseif nargin<3
        bias = 0;
    end

    % find the differences
    srcharray = srcharray-srchvalue;

    if bias == -1   % only choose values <= to the search value

        srcharray(srcharray>0) =inf;

    elseif bias == 1  % only choose values >= to the search value

        srcharray(srcharray<0) =inf;

    end

    % give the correct output
    if nargout==1 | nargout==0

        if all(isinf(srcharray(:)))
            r = [];
        else
            r = find(abs(srcharray)==min(abs(srcharray(:))));
        end 

    elseif nargout>1
        if all(isinf(srcharray(:)))
            r = [];c=[];
        else
            [r,c] = find(abs(srcharray)==min(abs(srcharray(:))));
        end

        if nargout==3
            V = srcharray(r,c)+srchvalue;
        end
    end
end