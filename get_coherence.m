% compute coherence using a weighted ratio (accouting for spatial distribution) of eigenvalues 
%
% INPUT
% mEnv  - time-series matrix (channel X timesample) 
%
% OUTPUT
% coh   - coherence as the weighted ratio of eigenvalues 
function coh = get_coherence(mEnv)
% zero mean
x = (mEnv-repmat(mean(mEnv,2),[1 size(mEnv,2)]))';
% correlation matrix
C = corr(x);

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