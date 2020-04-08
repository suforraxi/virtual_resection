% USAGE 
%    [S1,M] = transfun(S,g,L);
% FUNCTION:
%    transforms and interpolates S-> g(S) for values of S between [L(:,1) L(:,2)) bin walls
%    for values outside L, S1=0
% INPUT:
%   S - any 1d signal
%   g - n-point vector with function values
%   L - [n,2] array of bin walls
% OUTPUT:
%   S1 - the transformed signal
%   M - masc for the transformed points
% VERSION
%   Stiliyan 05.12.2007
% Copyright (C) 2020 Stiliyan Kalitzin
function [S1,M] = transfun_i(S,g,L)
n=size(L,1); 
S1=S*0;
M=logical(S*0); 
Lm=mean(L,2); 
k=1;
m = (S>=L(k,1))&(S<L(k,2));
S1(m)=(g(k)*(Lm(k+1)-S(m))+g(k+1)*(S(m)-Lm(k)))/(Lm(k+1)-Lm(k));
M(m)=true; 
for k=2:n-1
  m = (S>=L(k,1))&(S<L(k,2));
  S1(m)=g(k)*(S(m)-Lm(k-1)).*(S(m)-Lm(k+1))/(Lm(k)-Lm(k-1))/(Lm(k)-Lm(k+1))+ ...
        g(k+1)*(S(m)-Lm(k)).*(S(m)-Lm(k-1))/(Lm(k+1)-Lm(k))/(Lm(k+1)-Lm(k-1))+ ...
        g(k-1)*(S(m)-Lm(k)).*(S(m)-Lm(k+1))/(Lm(k-1)-Lm(k))/(Lm(k-1)-Lm(k+1));
  M(m)=true; 
end;
 k=n; 
 m = (S>=L(k,1))&(S<L(k,2));
 S1(m)=(g(k-1)*(Lm(k)-S(m))+g(k)*(S(m)-Lm(k-1)))/(Lm(k)-Lm(k-1));
 M(m)=true; 

