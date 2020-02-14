% USAGE 
%    [S1,M] = transfun(S,g,L);
% FUNCTION:
%    transforms S-> g(S) for values of S between [L(:,1) L(:,2)) bin walls
%    for values outside L, S1=0
% INPUT:
%   S - any 1d signal
%   g - n-point vector with function values
%   L - [n,2] array of bin walls
% OUTPUT:
%   S1 - the transformed signal
%   M - masc for the transformed points
% VERSION
%   Stiliyan 25.11.2004

function [S1,M] = transfun(S,g,L);
S1=S*0;
M=logical(S*0); 
n=size(L,1); 
for k=1:n
  m = (S>=L(k,1))&(S<L(k,2));
  S1(m)=g(k);
  M(m)=logical(1); 
end;
