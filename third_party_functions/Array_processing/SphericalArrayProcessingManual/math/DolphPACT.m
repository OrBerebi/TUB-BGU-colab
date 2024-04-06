function Z = DolphPACT(N);
%DolphPACT returns the matrix product Z = P*A*C*T related to the design of
% a Dolph-Chebyshev beam pattern.
% Z = DolphPACT(N);
% N is the order.
%
% Fundmentals of Spherical Array Processing
% Boaz Rafaely, 2017.

P = zeros(N+1,N+1);
A = P;
C = P;
T = P;
Z = P;

% Generate the P matrix
for i = (0:N);
    p = legendre_coefficients(i);
    p = p(end:-1:1);
    p = [p zeros(1,N-i)];
    P(i+1,:) = p;
end


% Create the A matrix
i = repmat((0:N)',1,N+1);
j = i';
A = 2./(i+j +1);
A(mod(i+j,2) ~=0) = 0;


% Generate the C matrix
C = factorial(j)./(factorial(i).*factorial(abs(j-i)).*(2.^j));
C(i>j) = 0;

% Create the T matrix
TT =chebyshev_coefficients(2*N);
TT = TT(1:2:end);
TT = TT(end:-1:1);
T(1:N+2:end) = TT;

Z = P*A*C*T;
