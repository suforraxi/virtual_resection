% Compute the envelope using Hilbert transform
%
% INPUT
% cfg     - configuration struct with the following fields
%           EnvType                       - type of computation for the envelope it could be 
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

%     Copyright (C) 2020 Matteo Demuru
% 
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <https://www.gnu.org/licenses/>.

function e_m = get_Envelope(cfg,m)
    
    e_m = zeros(size(m));
    switch cfg.EnvType
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