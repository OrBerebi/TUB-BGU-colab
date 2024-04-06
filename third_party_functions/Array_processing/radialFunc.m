function out = radialFunc (n,kr,sphereType,a_max,timeDep,ext)
% function out = radialFunc (n,kr,sphereType,a_max,timeDep,ext)
%
% Generate a radial function:
%   n           - order
%   kr          - k times sphere radius (r); measurement position is
%                 assumed to be on the sphere.
%
%   sphereType  - type of radial function:
%                   0 = Open sphere
%                   1 = Rigid sphere
%                   2 = Open sphere with Cardioid mics
%
%   [a_max]     - optional: soft-limiting factor (dB), {default = inf, i.e. no limiting}
%   [timeDep]   - optional: time dependence assumption, {default = true}
%                       true  = positive, i.e. exp(jwt)
%                       false = negative, i.e. exp(-jwt)
%
%   [ext]       - optional: process externally using SOFiA {default=false}
%                   if true, then the SOFiA toolbox must be installed and
%                   present in the MATLAB path(!)
%
%   Output      - 1 x length(kr) vector
%
% July 2014, Jonathan Sheaffer, Boaz Rafaely, Ben-Gurion University
% Part of the EARS Beamforming Toolbox
%


  if nargin < 6
      ext = false;
  end
  if nargin < 5
      timeDep = true;
  end
  if nargin < 4
      a_max = inf;
  end
  
  kr = kr + (kr==0).*1e-30;   % Replace zero entries in kr

  if ext == false     % Process internally
      inpart = 4*pi*1i^n;
      if sphereType==0          % open
          
          rF= inpart.*sBessel(n,kr);
          
      elseif sphereType==1      % rigid
          
          if timeDep, outpart=sBessel(n,kr)-(sBesseld(n,kr)./sHankeld(n,kr,2)).*sHankel(n,kr,2);
          else outpart=sBessel(n,kr)-(sBesseld(n,kr)./sHankeld(n,kr,1)).*sHankel(n,kr,1); end;
          rF = inpart.*outpart;
          
      elseif sphereType==2      % cardioid
          rF = inpart .* (sBessel(n,kr) - 1i*sBesseld(n,kr) );         
      end
      
      if a_max==inf         % no soft-limiting
          
          out = rF;
          
      else                  % soft-limit
          
          alp=10^(a_max/20);
          bnrat=abs(rF)./rF;
          bnpi=atan(pi./(2.*alp.*abs(rF)));
          out=1./((2*alp/pi)*bnrat.*bnpi);
          
      end
      
  else                % Process externally
    disp('External processing chosen. Ignoring timeDep variable.');  
    if sphereType==0, sofia_sph=0;
    elseif sphereType==1, sofia_sph=2;
    elseif sphereType==2, sofia_sph=1;
    end
    
    if a_max==inf
        [dn, beam] = sofia_mf(n, double(kr), sofia_sph);
    else
        [dn, beam] = sofia_mf(n, double(kr), sofia_sph, a_max);
    end
    
    out = 1./dn(end,:);
    
  end


end

% internal functions
function y = sBessel(n,x)
  y = sqrt(pi./(2*x)) .* besselj(n+0.5,x);
end

function y = sBesseld(n,x)
  y = (n./x) .* sBessel(n,x) - sBessel(n+1,x);
end

function y = sHankel(n,x,kind)  % kind is 1 or 2
  y = sqrt(pi./(2*x)) .* besselh(n+0.5,kind,x);
end

function y = sHankeld(n,x,kind)  % kind is 1 or 2
  y = (n./x) .* sHankel(n,x,kind) - sHankel(n+1,x,kind);
end

