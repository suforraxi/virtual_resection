% Compute the envelope using Hilbert transform
%
% INPUT
% m       - time-series matrix (channel X time)
%
% OUTPUT
% e_m     - matrix of envelopes (channel X time)

function e_m = get_Envelope(m)

    e_m = hilbert(m');
    e_m = abs(e_m)';
    
%     e_m = zeros(size(m));
%     for i = 1 : size(m,1)
%     
%         s        = huang(m(i,:),1,2);
%         e_m(i,:) = s(2,:); 
%     end