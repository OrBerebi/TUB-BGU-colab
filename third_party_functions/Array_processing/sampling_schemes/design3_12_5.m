function [a,th,ph] = design3_12_5;

% function [theta phi] = design3_12_5;
% a - weights (constant)
% th- theta for all points
% ph - phi for all points
%
% Spherical design, from Hardin and Sloane
% http://www.research.att.com/~njas/sphdesigns/dim3/des.3.12.5.txt
% Design for 5 order polynomial, i.e. n<=2 in our case
% uses 12 samples on the sphere
% This uses regular isocahedron configuration

j=1:12;

xyz=[
0.850650808352E+00
0			
-0.525731112119E+00
0.525731112119E+00
-0.850650808352E+00
0.000000000000E+00
0		
-0.525731112119E+00
0.850650808352E+00
0.850650808352E+00
0			
0.525731112119E+00
-0.525731112119E+00
-0.850650808352E+00
0
0		
0.525731112119E+00
-0.850650808352E+00
-0.850650808352E+00
0			
-0.525731112119E+00
-0.525731112119E+00
0.850650808352E+00
0
0		
0.525731112119E+00
0.850650808352E+00
-0.850650808352E+00
0			
0.525731112119E+00
0.525731112119E+00
0.850650808352E+00
0
0		
-0.525731112119E+00
-0.850650808352E+00
]';

x=xyz(3*(j-1)+1);
y=xyz(3*(j-1)+2);
z=xyz(3*(j-1)+3);

%r=sqrt(x^2+y^2+z^2);
th=acos(z);
ph=atan2(y,x);

a=ones(size(th))*(4*pi/12);