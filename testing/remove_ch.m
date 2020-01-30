% Apply orthogonalization to  data
%
% INPUT
% d      - fieldtrip data structure
% idx_ch - index of the channel to be removed from the data 
%
% OUTPUT
% d      - fieldtrip data structure after the orthogonalization

function d = remove_ch(d,idx_ch)

    nTrial      = numel(d.trial); 
    idx         = true(length(d.label),1);
    idx(idx_ch) = false;
    
    for i = 1 : nTrial 
        m          = d.trial{i};
        m          = get_ortho_matrix(m,idx_ch);
        d.trial{i} = m;
        d.label    = d.label(idx);   
        
    end