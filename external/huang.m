% USAGE
%  [S1,N,D,F,P] = huang(S,n1,n2,method);
% FUNCTION:
%  gives Huang decomposition of a 1D signal
% INPUT:
%  S - input 1D signal
%  n1: if n1>=1 number of iterations for each S1 cmponent, 
%      if 0<n1<1, then std(W1)/std(R1)<n1 is reached
%      if -1<n1<0, iteration stops when |sum(extrema)|/sum(|extrema|)<-n1
%      if n1<=-1, iteration stops when the number of signal maxima that are
%                               <0 plus minima that are >0 is less than -n1
%      default n1=-1; (the smart choice)  
%  n2: number of decomposition steps = size(S1,2), 
%      if 0<n2<1, then std(R1)/std(S)<n2 is reached
%      if n2<=0, then #of extema<=-n2 of the last component is reached
%      default n2=-1; (the smart choice) 
%  method: string with interpolation method, 
%      default is 'interflat'
% OUTPUT:
%  S1 - 2D array with Huang components S1(component, sample); 
%       last component is the residual
%  N -  the number of extrema per component
%  D- 2D array with envelopes (max-min)/2 of the huang components (component, sample); 
%  F- 2D array with the instanteneous frequencies of the huang components (component, sample); 
%  F- 2D array with the instanteneous phases of the huang components (component, sample); 
% USES
%  get_extrema
%  interflat
% VERSION:
%  Stiliyan & George, 03.08.2017
% Copyright (C) 2020 Stiliyan Kalitzin
function [S1,N,D,F,P] = huang(S,N1,N2,method)
if nargin<2
    N1=[]; 
end
if nargin<3
    N2=[]; 
end
if nargin<4
    method=[]; 
end
ph=nargout>3;
if isempty(method)
    method='interflat';
end
if isempty(N1)
    N1=1; 
end
if isempty(N2)
    N2=1; 
end
stop=false;
k=0; 
R1=S; 
while ~stop
    k=k+1; 
    if ph
   [S1(k,:),R1,N(k),huj,D(k,:),F(k,:),P(k,:)]=huang_1(R1,N1,method);
    else
    [S1(k,:),R1,N(k),huj,D(k,:)]=huang_1(R1,N1,method);
    end      
    if N2>=1 
         stop=k>=N2; 
    elseif N2>0
         stop=std(S1(k,:))/std(S)<=N2; 
    else
        stop=N(k)<=-N2;
    end
%      if (k>2)&(N(k-1)<=N(k))
%          stop=true;
%          N=N(1:end-1); 
%          S1=S1(1:end-1,:); 
%      end;
end
S1=[S1;R1]; 
n=size(S1,1); 
[JJ,UU]=get_extrema(R1); 
N=[N numel(JJ)]; 
if ~ph
    return;
end
if numel(JJ)==0
    F(n,:)=0;
    P(n,:)=0;
else
    for e=1:numel(JJ)-1
        F(n,JJ(e):(JJ(e+1)-1))=1/(JJ(e+1)-JJ(e));
    end
    F(n,1:JJ(1))=F(JJ(1));
    F(n,JJ(end):end)=F(JJ(end));
    P(n,:)=S*0;
    for e=1:numel(JJ)-1
        P(n,JJ(e):(JJ(e+1)-1))=-pi*(1+sign(UU(e)))/2+pi*(0:(JJ(e+1)-JJ(e)-1))/(JJ(e+1)-JJ(e));
    end
end


function [S1,R1,u,k,D,F,P] = huang_1(S,N,method) 
pf = nargout>5; 
k=0; 
stop=false; 
S1=S;
R1=S1*0; 
W1=S;
D=S*0; 
F=S*0; 
P=S*0; 
if std(S)==0
    u=0; 
    return; 
end 
while ~stop
     k=k+1; 
     [JJ,UU,AA]=get_extrema(S1);
     if isempty(JJ)
         break;
     end
     % stop criterium
     if N>=1
         stop=k>N;
     elseif N>0
         stop=std(W1)/std(S)<=N;
     elseif N<=-1
         q=sign(UU).*sign(AA);
         stop=(sum(q==1)<-N);
     else
         q=abs(sum(AA))/sum(abs(AA)); 
         stop=q<-N; 
     end
     J=[1 JJ length(S1)];
     A=[S1(1) AA S1(end)];
     if numel(UU)>0
         U=[-UU(1) UU -UU(end)];
     else
         U=[-2 2]*sign(S1(1)-S1(end));
     end
     if stop
         break;
     end
     
     % interpolation of the maxima
     j1=J(U==-2);
     A1=A(U==-2);
     if j1(1)>1 
         j1=[1 j1];
         A1=[A1(1) A1]; 
     end
     if j1(end)<length(S1)
         j1=[j1 length(S1)]; 
         A1=[A1 A1(end)]; 
     end
     if strcmp(method,'interflat') 
        M=interflat(j1, A1,[1:length(S1)]); 
     else
        M=interp1(j1, A1,[1:length(S1)],method,'extrap');      
     end
   
% interpolation of the minima
    j2=J(U==2); 
    A2=A(U==2); 
      if j2(1)>1 
         j2=[1 j2];
         A2=[A2(1) A2]; 
      end
     if j2(end)<length(S1)
         j2=[j2 length(S1)]; 
         A2=[A2 A2(end)]; 
     end
     if strcmp(method,'interflat') 
        m=interflat(j2, A2,[1:length(S1)]); 
     else
         m=interp1(j2,A2,[1:length(S1)],method,'extrap');
     end
     % updating the data
       W1=(M+m)/2;
       S1=S1-W1; 
       R1=R1+W1;
       if k==1
          D=(M-m)/2;
       end
end
[JJ,UU]= get_extrema(S1);
u=numel(JJ); 
if ~pf
    return;
end
% instan frequencies & phases
if u>0
    for e=1:numel(JJ)-1
        F(JJ(e):(JJ(e+1)-1))=1/(JJ(e+1)-JJ(e))/2;
    end
    F(1:JJ(1))=F(JJ(1));
    F(JJ(end):end)=F(JJ(end));
    for e=1:numel(JJ)-1
        P(JJ(e):(JJ(e+1)-1))=-pi*(1+sign(UU(e)))/2+pi*[0:(JJ(e+1)-JJ(e)-1)]/(JJ(e+1)-JJ(e));
    end
end
return;





