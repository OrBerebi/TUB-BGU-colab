function M=getRotationMat(psi,theta,phi)
%calculating rotatin matrix
%The rotation axis sequences is Z-X-Z and the rotation direction is counterclockwise
% with ?, ?, ? as the respective angles. The rotation matrix is M = R1(?)*R2(?)*R3(?)

R1=eye(3);
R1(1:2,1:2)=[cos(phi) -1*sin(phi) ;sin(phi) cos(phi) ];

R2=eye(3);
R2(2:3,2:3)=[cos(theta) -1*sin(theta) ;sin(theta) cos(theta) ];

R3=eye(3);
R3(1:2,1:2)=[cos(psi) -1*sin(psi) ;sin(psi) cos(psi) ];

M=R1*R2*R3;
