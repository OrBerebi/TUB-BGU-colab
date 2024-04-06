function out=meanRespN(N, kr, timeDep, useSOFiA, useTaperWin, winVec)
% function out=meanRespN(N, kr)
%
% Calculate averaged magnitude pressure at left and right ears
% due to diffraction around a rigid sphere, for order N

if ~exist('useTaperWin','var')
    useTaperWin = false;
    
end
bn=radialMatrix (N,kr,1,inf,timeDep,useSOFiA);
Yl = shMatrix(N,pi/2,pi/2);

if useTaperWin
    tmp = (winVec.*abs(bn).^2)*(abs(Yl).^2);
else
    tmp=(abs(bn).^2)*(abs(Yl).^2);
end

out=sqrt(tmp/(4*pi));


end