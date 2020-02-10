% Compute the envelope using Hilbert transform
%
% INPUT
% cfg     - configuration struct with the following fields
%           type                       - type of computation for the envelope it could be 
%
%                                       'filt_hilbert_env' assuming the time-series is band filtered it computes the envelope using Hilbert
%
%                                       'huang_env'        compute the empirical mode decomposition (see huang.m) and then select (cfg.seleted_comp) 
%                                                          the component for which the envolope will be computed using Hilbert 
%           nHComp                     - should be set when using huang_env 
%                                        number of huang component to consider (iterative subdivision of the signals)  
%           selected_comp              - should be set when using huang_env 
%                                        index of the component for which the envelope will be computed  
%           
% m       - time-series matrix (channel X time)
%
%
% OUTPUT
% e_m     - matrix of envelopes (channel X time)

function e_m = get_Envelope(cfg,m)
    
    e_m = zeros(size(m));
    switch cfg.type
        case 'filt_hilbert_env'    
              e_m = hilbert(m');
              e_m = abs(e_m)';
        case 'huang_env'
            
            nHComp        = cfg.nHComp; 
            selected_comp = cfg.selected_comp; 
            
            for i = 1 : size(m,1)

                s        = huang(m(i,:),1,nHComp);
                comp     = s(selected_comp,:);
                comp     = abs(hilbert(comp'))';
                e_m(i,:) = comp; 

            end
    end