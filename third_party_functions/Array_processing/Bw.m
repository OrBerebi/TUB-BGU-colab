function y = Bw(n,kr,ka,sphere);

% Equation according to Williams 1999
% with k denoting arrival direction exp(i[wt+kr])
%
% y = Bw(n,kr,ka,sphere);
% Bn(kr,ka), order n
% sphere radius a
% sphere=0 open  sphere
% sphere=1 rigid sphere
% sphere=2 cardiod mic open sphere
% the term (-j)^n is used so that (theta,phi) are direction of arrival (not propogation) 

j=sqrt(-1);

if sphere==0, %open
   y = 4 * pi * j^n * besseljs(n,kr);

elseif sphere==1, %rigid
   y = 4 * pi * j^n * ( besseljs(n,kr) - ( besseljsd(n,ka)./besselhsd(n,ka) ) .* besselhs(n,kr) );
   
elseif sphere==2, %cardiod
   y = 4 * pi * j^n * ( besseljs(n,kr) - j*besseljsd(n,kr) );
   
else,
    y=0;
    
end
