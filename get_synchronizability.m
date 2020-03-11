% compute global coherence using a weighted ratio (accounting for spatial distribution) of eigenvalues 
%
% INPUT
%
% C     - symmetric functional connectivity matrix 
%
% OUTPUT
% 
% syn   - 'synchronizability' as the ratio of weighted (using eigenvectors) eigenvalues
%              
% V     - eigenvectors matrix
% E     - eigenvalues  matrix

function [syn, V, E] = get_synchronizability(C)


% eigen-value decopmposition
[V,E] = eig(C); 

% weighted ratio
%f    = E*(mean(V.^2)./max(V.^2))';
%syn = max(f)/sum(f);


% kini 2019

ord_E = sort(diag(E),'ascend');

% second smallest eigenvalue over max eigenvalue
syn   = ord_E(2)/ord_E(end);