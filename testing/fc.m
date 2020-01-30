
% compute functional connectivity
%
% INPUT 
% x     - data matrix (channel X timesample)
%
% OUTPUT 
% 
% C     - functional connectivity matrix (correlation) 
%
function C = fc(x)

    C = corr(x');

