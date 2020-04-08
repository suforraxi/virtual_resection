% USAGE:
%    R2 = h2_remove(S2,S1,n)
% FUNCTION:
%    removes the S1 nonlinear component from S2
% INPUT:
%    S1 - driving signal
%    S2 - target signal
%    n - bin size, the same as in h2
% OUTPUT:
%    R2 - the cleaned S2 signal
% VERSION:
%    Stiliyan, 05.12.2007
% Copyright (C) 2020 Stiliyan Kalitzin

function R2 = h2_remove(S2,S1,n)
if nargin<3 
    n=[]; 
end; 
[h,g,N,D,L] = h2(S2,S1,n);
m=(isfinite(g));
L=L(m,:);
g=g(m); 
if size(L,1) > 1 
    R2= S2-transfun_i(S1, g ,L);
else
    R2= S2-transfun(S1, g ,L);
end;
