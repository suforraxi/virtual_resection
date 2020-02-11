
% compute functional connectivity
%
% INPUT 
% x       - data matrix (channel X timesample)
% fc_type - string one of the possible
%          'corr'
%          'coh'
% OUTPUT 
% 
% C     - functional connectivity matrix (correlation) 
%
function C = fc(x,fc_type)

nch = size(x,1);
C   = zeros(nch,nch);

switch fc_type
    case 'corr'
        C              = corr(x');
        %C(1:nch+1:end) = 1;
        %C              = diag(sum(C,2))-C;
    case 'coh'
        % a matrix (channels X time)
    
     
        h          = hilbert(x'); 

        for i = 1 : nch
            for j = 1 : nch
                if i<j
                   
                   Am     = abs(h(:,i)).*abs(h(:,j));
                   phDiff = exp(1i * angle( ( h(:,i).*conj(h(:,j)) )));
                   C(i,j) = abs(mean(Am.*phDiff)/mean(Am));
                end 
            end 
        end
        
       
        C              = C + C';
        C(1:nch+1:end) = 1;
        %C              = diag(sum(C,2))-C;

    case 'xcorr'
        
        C = xcorr(x',1024,'normalized');
        C = max(abs(C))';
        C = reshape(C,size(x,1),size(x,1));
         

end


