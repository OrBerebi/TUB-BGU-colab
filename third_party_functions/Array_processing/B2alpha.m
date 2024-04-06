function BB = B2alpha (N,kr,sphere, alpha)
% same as B2 only with soft limiting (alpha)

    if sphere==0, sofia_sph=0;
    elseif sphere==1, sofia_sph=2;
    elseif sphere==2, sofia_sph=1;
    end
    
    if nargin<4
        [dn, beam] = sofia_mf(N, double(kr), sofia_sph);
    else
      [dn, beam] = sofia_mf(N, double(kr), sofia_sph, alpha);
    end
    dn=1./dn.';     % turn d_n into b_n
    count=2;
    Dn(:,1)=dn(:,1);
    for nn=1:N
      for mm=-nn:nn
        Dn(:,count) = dn(:,nn+1);  
        count=count+1;
      end
    end
    
    BB=Dn;
    
end