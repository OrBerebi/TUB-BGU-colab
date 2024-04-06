function BB = BnMat(N,kr,ka,sphere);
%BnMat return the radial functions in matrix form.
% y = BnMat(n,kr,ka,sphere);
% n is the order.
% k is the wave number.
% a is the radius of the open/rigid sphere.
% r is the radius of the microphone positions.
% sphere options:
%  0 open sphere with pressure micropnones.
%  1 rigid sphere with pressure microphones.
%  2 open sphere with Cardiod microphones.
%
% Fundmentals of Spherical Array Processing
% Boaz Rafaely, 2017.

% complex constant
j=sqrt(-1);
BB=[];

for n=0:N,
    B = Bn(n,kr,ka,sphere);
    BB = [BB; repmat(B,2*n+1,1)];
end;

BB=BB.';