% S1 =  h2_partialize(S,IX,n);
function S1 = h2_partialize(S,IX,n)
if nargin<3
    n=0.1;
end
J=(1:size(S,1)); 
for k=1:numel(IX)
    c=IX(k); 
    J=setdiff(J,c); 
    for k1=1:numel(J)
        c1=J(k1); 
        S(c1,:)=h2_remove(S(c1,:),S(c,:),n); 
    end
end
S1=S(J,:); 