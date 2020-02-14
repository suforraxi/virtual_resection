% playing with h2

function playing_h2()

close all

t  = 1 : 1000;
fr = [11 16 25 ];

x = zeros(numel(fr),numel(t));

for i = 1 : numel(fr)
    x(i,:) = sin(2*pi*fr(i)*t./length(t));
end

y = zeros(numel(fr),numel(t));
for i = 1 : 2
    y(i,:) = x(i,:) + x(end,:);
end

y(end,:) = x(end,:);



C = compute_h2(y);

figure
imagesc(C,[0 1])
colormap('jet')
colorbar;


s1 = y(1,:) ;   
s2 = y(2,:) ;

[m,u] = h2_local(s1,s2,40) ;

function m = compute_h2(x)
    
    m = zeros(size(x,1));
    
    for i = 1 : size(x,1)
    
        for j = 1 : size(x,1)
            
            m(i,j) = h2_m(x(i,:),x(j,:));
        
        end
    end


function [C,u]=h2_local(S1,S2,n) 
if n>0 
    S2=floor(n*(S2-min(S2))/(max(S2)-min(S2)))+1;
    S2(S2>n)=n;
else 
    n=max(S2); 
end
D=zeros(1,n); 
M=zeros(1,n); 
N=zeros(1,n); 
for k=1:length(S2) 
  D(S2(k))=D(S2(k))+S1(k).^2;
  M(S2(k))=M(S2(k))+S1(k);
  N(S2(k))=N(S2(k))+1; 
end
D(N>0)=D(N>0)-M(N>0).^2./N(N>0); 
D=sum(D); 
C=1-D/sum(sum(S1.^2)); 
u=max(N)/sum(N)-1/n; 