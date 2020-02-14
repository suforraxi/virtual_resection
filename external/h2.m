% USAGE:
%   [C,g,N,D,L]=h2(S1,S2,n);
% FUNCTION:
%   gives the non-linear assymetric association h^2 index b/n signals S1 and S2, 
%   S2->S1 = var(S1|S2)/var(S1)
% INPUT:
%    S1, S2: 1D signals
%    n : number of bins for S2, default [100]
%        if n<=1, the bin size is that fraction of the std(S2)
%        if n<0 the S2 interval is devided into |n| bins of equal # of points
%        if n is a vector it defines the bin's borders
% OUTPUT:
%  C - h^2 association index = (total.variance-explained.variance)/total.variance
%  g - the regression function = mean(S1|S2), runs from S2min to S2max
%  N - nimber of points in the k-th bin 
%  D - deviation in the k-th bin std(S1|S2)
%  L - 2D array with bin borders of S1; [L(:,1) L(:,2)]
% VERSION:
%  Stiliyan, 05.08.05

function [C,g,N,D,L]=h2(S1,S2,n);
S1=reshape(S1,[1 numel(S1)]); 
S2=reshape(S2,[1 numel(S2)]); 
C=0;
if nargin <2  return; end;
if nargin<3 n=[]; end; 
if isempty(n) n=100; end;
m=min(length(S1),length(S2));
if m==0 return; end;
S1=S1(1:m);
S2=S2(1:m);
m=isnan(S1)|isnan(S2)|isinf(S1)|isinf(S2);
S1=S1(~m); 
S2=S2(~m); 
if length(S1)<1 retun; end; 
if n<=1 n=round((max(S2)-min(S2))/(n*std(S2))); end; 
A1=mean(S1); 
A=sum((S1-A1).^2);
if A==0  return; end;
d=max(S2)-min(S2);
if d==0 return; end;
if length(n)==1
    if n>0 
        L=linspace(min(S2),max(S2),n+1);
    else 
        n=-n;
       P=linspace(0,100,n+1);
       for k=1:length(P) L(k)=prctile(S2,P(k)); end;
   end;
else
    L=n; 
    n=length(L)-1; 
end;
N=zeros(1,n); 
g=zeros(1,n); 
D=zeros(1,n); 
for k=1:n
   M=S1((S2>=L(k))&(S2<L(k+1)));
   N(k)=length(M);
   if N(k)>0 
       g(k)=mean(M); 
       D(k)=sum((M-g(k)).^2); 
   end;
end;
B=sum(D(N>0)); 
 C=(A-B)/A;
D(N>0)=sqrt(D(N>0)./N(N>0)); 
L=[L(1:end-1)' L(2:end)'];