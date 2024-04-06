function T = chebyshev_coefficients(N);
%CHEBYSHEV_COEFFICIENTS Returns the coefficients of Chebyshev polynomials. 
% T = chebyshev_coefficients(N);
% N is the order.
% It uses the recursion Tn+1(x) = 2 x Tn(x) - Tn-1(x).
%
% Fundmentals of Spherical Array Processing
% Boaz Rafaely, 2017.

if N==0,
    T=1;
    
elseif N==1,
    T=[1 0];
    
else  
  T0=1;
  T1=[1 0];
  for n=2:N,
      T2=2*[T1 0]-[0 0 T0];
      T0=T1;
      T1=T2;
  end;
  T=T1;
end;
  