% compute the Laplacian Matrix from a functional connectivity matrix
% 
% INPUT
%
% C     - symmetric functional connectivity matrix
%
% OUTPUT
%       
% L     - laplacian matrix (degree matrix - functional connectivity matrix)

function L = get_laplacian(C)

  L = diag(sum(C,2)) - C;