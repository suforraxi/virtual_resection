
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


