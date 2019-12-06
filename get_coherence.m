% compute coherence using a weighted ratio (accouting for spatial distribution) of eigenvalues 
%
% INPUT
%
% C     - symmetric functional connectivity matrix 
%
% OUTPUT
% 
% coh   - coherence as the ratio of weighted (using eigenvectors) eigenvalues 
% V     - eigenvectors matrix
% E     - eigenvalues  matrix
function [coh, V, E] = get_coherence(C)

% eigen-value decopmposition
[V,E] = eig(C); 


nEig = size(V,1);
%f    = zeros(1,nEig); 
f    = E*(mean(V.^2)./max(V.^2))';


%  for i = 1 : nEig
%  
%     f(i) = E(i,i) * (mean(V(:,i).^2)/max(V(:,i).^2)); 
%  end
%  
coh = max(f)/sum(f);

