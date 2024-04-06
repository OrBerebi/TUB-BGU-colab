function [a,th,ph]=uniform_sampling_extended(N)

% [a,th,ph]=uniform_sampling(N);
%
% a - weights
% th - theta angles for all points
% ph - phi angles for all points
% Total no. of points is length(a)=length(th)=length(ph)
%
% Nearly uniform sampling for order N
%
% Boaz Rafaely 23 June 2004
%
% new designs added by David Alon 09 Nov 2010
%new designs of order 7 and 9
%based on "McLaren’s Improved Snub Cube and Other New Spherical Designs in Three Dimensions"
%  R. H. Hardin and N. J. A. Sloane
% all are orthogonal designs





if N==2 % Total 12 samples
    [a,th,ph]=design3_12_5;

elseif N==3 % Total 32 samples
    [a,th,ph]=design3_32_7;
    
elseif N==4 % Total 36 samples
    [a,th,ph]=design3_36_8;
    
elseif N==5 % Total 73 samples
%     [a,th,ph]=design3_64_10;
    [a,th,ph]=design3_73_10;

elseif N==6 % Total 84 samples
    [a,th,ph]=design3_84_12;

elseif N==7 % Total 108 samples
    [a,th,ph]=design3_108_14;
    
elseif N==8 % Total 144 samples
    [a,th,ph]=design3_144_16;

elseif N==9 % Total 180 samples
    [a,th,ph]=design3_180_18;
    
elseif N==10 % Total 240 samples
    [a,th,ph]=design3_240_21;
    
else
    fprintf('\n Not available for given N');
    return;
end;


