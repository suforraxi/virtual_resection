% Apply non linear partialization to data
%
% INPUT
% d      - fieldtrip data structure
% idx_ch - indexes of the channels to be removed from the data 
%
% OUTPUT
% d      - fieldtrip data structure after the non linear partialization

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

function d = removeNoNLinear_ch(d,idx_ch)

    nTrial      = numel(d.trial); 
    idx         = true(length(d.label),1);
    idx(idx_ch) = false;
    
    for i = 1 : nTrial 
        m          = d.trial{i};
        m          = h2_partialize(m,find(idx_ch));
        d.trial{i} = m;
          
    end
     d.label    = d.label(idx); 