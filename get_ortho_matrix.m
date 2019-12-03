% orthogonalization of signals using Gram–Schmidt
% 
%
% INPUT 
% m       time-series matrix (channel X time) size (M,N)
% idx_ch  index of the channel to use for the orthogonalization   
%
% OUTPUT
% o_m     orthogonalized time-series (channel X time) size (M-1,N)


function o_m = get_ortho_matrix(m,idx_ch)

    idx2keep         = true(1,size(m,1));
    idx2keep(idx_ch) = 0;
    u                = m(idx_ch,:);
    o_m              = m(idx2keep,:);
    nCh              = size(o_m,1);
    
    
    X   = o_m;
    U   = repmat(u,[nCh 1]);
    o_m = o_m - ( ( diag(X*U')/(u*u') ).*U );
    
    
%     for i = 1 : nCh
%         x        = o_m(i,:);
%         o_m(i,:) = get_ortho(x,u);
%     end
%     
    
    
     
    
function x_ortho = get_ortho(x,y)
    
    x_ortho = x - ((x*y')/(y*y'))*y;