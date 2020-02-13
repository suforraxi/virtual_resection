% Apply orthogonalization to  data
%
% INPUT
% d      - fieldtrip data structure
% idx_ch - index of the channel to be removed from the data 
%
% OUTPUT
% d      - fieldtrip data structure after the non linear orthogonalization

function d = removeNoNLinear_ch(d,idx_ch)

    nTrial      = numel(d.trial); 
    idx         = true(length(d.label),1);
    idx(idx_ch) = false;
    
    for i = 1 : nTrial 
        m          = d.trial{i};
        m          = h2_partialize(m,idx_ch);
        d.trial{i} = m;
          
    end
     d.label    = d.label(idx); 