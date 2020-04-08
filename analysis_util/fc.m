
% compute functional connectivity
%
% INPUT 
% x       - data matrix (channel X timesample)
% fc_type - string one of the possible measures
%          'corr'  - Pearson correlation (zero lag)
%          'xcorr' - cross-correlation  (choose max value across delays)
%          'coh'   - coherence (see spectral coherence from Analyzing Neural Time Series Data Mike X Cohen, MIT press)
%          'h2'    - non-linear correlation (see Kalitzin 2007 IEEE Trans Biom Eng)
%          'pli'   - phase lag index (see Stam 2007 Hum Brain Mapp)
% OUTPUT 
% 
% C     - functional connectivity matrix 
%
function C = fc(x,fc_type)

nch = size(x,1);
C   = zeros(nch,nch);

switch fc_type
    case 'corr'
        C              = corr(x');
      
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
       

    case 'xcorr'
        
        C = xcorr(x',1024,'normalized');
        C = max(abs(C))';
        C = reshape(C,size(x,1),size(x,1));
    
    case 'h2' % linear and non-linear
            
        for i = 1 : nch
            for j = 1 : nch
               
        
                  C(i,j) = h2_m(x(i,:),x(j,:));
            
            end 
        end
    
    case 'pli' 
    
        h          = hilbert(x'); 

        for i = 1 : nch
            for j = 1 : nch
                if i<j

                    C(i,j)  = abs( mean(sign( (angle( h(:,i).*conj(h(:,j)) ) ) ) ) ) ;

                end 
            end 
        end

        C              = C + C';
        C(1:nch+1:end) = 1;
    
end


