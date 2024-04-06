function [ degree ] = Order_to_degree_lebedev( N )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

%Allowed values are (degree -> order):
% order avaliable: { 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,17,20,
%                   23,26,29,32,35,38,41,44,47,50,53,56,59,62,65};

degrees_avail=[6, 14, 26, 38, 50, 74, 86, 110, 146, 170, 194, 230,... 
266, 302, 350, 434, 590, 770, 974, 1202, 1454, 1730, 2030, 2354, 2702,... 
3074, 3470, 3890, 4334, 4802, 5294, 5810];

orders_avail=floor(sqrt(degrees_avail/1.3)-1);

if isempty(find(orders_avail==N,1))
   error('Invalid quadrature degree');
end

degree = degrees_avail(find(orders_avail==N,1));

end

