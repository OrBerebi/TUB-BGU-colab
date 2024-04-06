function out = radialMatrix (N,kr,sphereType,a_max,timeDep,ext)
% function out = radialMatrix (N,kr,sphereType,a_max,timeDep,ext)
%
% Generate a matrix of radial functions:
%   N           - Maximum order
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
%   Output      - length(kr) x (N+1)^2 matrix
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


    if ext == false         % process internally

        bMat = [];
        for n=0:N
            bn = radialFunc(n,kr,sphereType,a_max,timeDep,false);
            bMat = [bMat; repmat(bn,2*n+1,1)];
        end

        out = bMat.';

    else                    % process externally
        
        disp('External processing chosen. Ignoring timeDep variable.');  
        if sphereType==0, sofia_sph=0;
        elseif sphereType==1, sofia_sph=2;
        elseif sphereType==2, sofia_sph=1;
        end
        
        if a_max==inf
            [dn, beam] = sofia_mf(N, double(kr), sofia_sph);
            %dn = sofia_rfi(dn);
        else
            [dn, beam] = sofia_mf(N, double(kr), sofia_sph, a_max);
            %dn = sofia_rfi(dn);
        end
        
        bn=1./dn.';     % turn d_n into b_n
        count=2;
        Bn(:,1)=bn(:,1);
        for nn=1:N      % reorder matrix according to our convention
          for mm=-nn:nn
            Bn(:,count) = bn(:,nn+1);  
            count=count+1;
          end
        end

        out=Bn;
        
    end


end