function y = Bn(n,kr,ka,sphere);
%BN returns the radial functions for a plane wave around a sphere.
% y = Bn(n,kr,ka,sphere);
% n is the order.
% k is the wave number.
% a is the radius of the open/rigid sphere.
% r is the radius of the microphone positions.
% sphere options:
%  0 open sphere with pressure micropnones.
%  1 rigid sphere with pressure microphones.
%  2 open sphere with Cardiod microphones.
%
% A model for the sound field is used that employs wave arrival directions
% and spherical Hankel functions of the second kind. For details see:
% Tourbabib and Rafaely, ACTA ACUSTICA UNITED WITH ACUSTICA, 101
% (2015),470-473.
%
% Fundmentals of Spherical Array Processing
% Boaz Rafaely, 2017.

j=sqrt(-1);

if sphere==0, %open
   y = 4 * pi * j^n * besseljs(n,kr);
                                                                                                                                               
elseif sphere==1, %rigid
   y = 4 * pi * j^n * ( besseljs(n,kr) - ( besseljsd(n,ka)./conj(besselhsd(n,ka)) ) .* conj(besselhs(n,kr)) );
   
elseif sphere==2, %cardiod
   y = 4 * pi * j^n * ( besseljs(n,kr) - j*besseljsd(n,kr) );
   
else,
    y=0;
    
end
