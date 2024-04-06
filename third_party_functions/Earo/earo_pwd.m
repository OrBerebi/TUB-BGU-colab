function [varargout] = earo_pwd(eobj, N, a_max, fVec)
% function [anm, [A]] = earo_pwd(eobj, N, [a_max])
%  
% Perform Plane-Wave Decomposition on an EARO object
%
%   eobj        - Input EARO object
%   N           - Decomposition order
%   [a_max]     - (optional) radial filter soft-limiting factor in dB, 
%                               {default = inf = no limiting}
%
%   [fVec]      - Frequency vector (must be supplied if eobj is in 'FREQ'
%                   domain!)
%
%   Output:
%   anm         - length(f) x (N+1)^2 matrix of plane wave density in the
%                   SH domain.
%   [A]         - (optional) length(f) x length(az/el) matrix of plane wave
%                   density coefficients in the space domain; computed over
%                   the azimuth/elevation grid specified in eobj.micGrid
%
% This function operates over the eobj.micGrid dimension.
%
% August 2014, Jonathan Sheaffer and David Alon, Ben-Gurion University
% Part of the EARO Beamforming toolbox
%

  useSOFiA = false;                     % change to true if you wish to generate the radial functions using SOFiA toolbox
                                        % instead of EARO's internal func.
  timeDep = true;                       % Time dependence assumption (true = positive, i.e. e^jwt); not available with SOFiA(!)

  if nargin<3, a_max=inf; end;          % Default = no soft-limiting
  

  % Transform to frequency domain, if needed
  if strcmp(eobj.dataDomain{1},'TIME')
      [eobj,fVec]=eobj.toFreq();
      if ~eobj.shutUp, disp('earo_pwd: Transforming to frequency domain.'); end;
  else
      if ~exist('fVec'), error('ERROR! eobj is in FREQ domain, but you did not supply a frequency vector!'); end;
  end

  % Transform to SH domain, if needed
  if strcmp(eobj.dataDomain{2},'SPACE')
      eobj=eobj.toSH(N,'MIC');
  elseif strcmp(eobj.dataDomain{2},'SH')
      if ~(eobj.orderN==N)  % if N and orderN are not the same, then fix it
          eobj=eobj.toSpace('MIC');
          eobj=eobj.toSH(N,'MIC');
      end
  end
  
  % Compute speed of sound
  if ~isempty(eobj.avgAirTemp) && isnumeric(eobj.avgAirTemp)
    c = 331.3 * sqrt(1+(eobj.avgAirTemp/273.15));
  else
      c = 343.5; % default value;
  end
  
  % Construct kr vector
  fVec=double(fVec);
  kr=double(fVec*2*pi*eobj.micGrid.r(1)/c);
  
  % Generate radial functions
  if eobj.scatterer, sphereType=1; else sphereType=0; end;
  bn = radialMatrix(N,kr,sphereType,a_max,timeDep,useSOFiA);
  
  % Perform PWD
  anm = double(squeeze(eobj.data(1,:,:)))./bn;
  
  % Output
  varargout{1} = anm;
  
  if nargout>1          % transform to space domain as well
      Y = shMatrix(N,eobj.micGrid.elevation,eobj.micGrid.azimuth); 
      varargout{2} = anm*Y;
  end
end

