function BB = B2(N,kr,ka,sphere)

% N maximum order for n
% kr vector for product of wave-number and radius

% complex constant
j=sqrt(-1);

BB=[];

for n=0:N,
    B = Bw(n,kr,ka,sphere);
    BB = [BB; repmat(B,2*n+1,1)];
end;

BB=BB.';