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
    