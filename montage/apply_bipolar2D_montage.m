% apply bipolar montage for grid and strips if present
% Assumption grid are 5x4 and strip 1x8 (max length)
% Grid montage is computed separately from strip montage and they are
% merged at the end (see merge_dateset)
%
% INPUT
% data    - fieldtrip data structure (see http://www.fieldtriptoolbox.org/)
%          
%             data.trial
%             data.time
%             data.fsample
%             data.label
%             data.sampleinfo
% OUTPUT       
% outdata - fieldtrip data structure with the data transformed according to the montage
%
function outdata = apply_bipolar2D_montage(data)
 
ntrial = numel(data.trial);

gridX = data;
gridY = data;
strip = data;
for i = 1: ntrial
    
    [gridX.label,gridX.trial{i}] = create_bipolar_montage_5by4X(gridX.label,gridX.trial{i});
  
    if(isempty(gridX.trial{i}))
        gridX = [];
        break
    end
    
end
for i = 1: ntrial
    
    [gridY.label,gridY.trial{i}] = create_bipolar_montage_5by4Y(gridY.label,gridY.trial{i});
     
    if(isempty(gridY.trial{i}))
        gridX  = [];
        break
    end
    
end
for i = 1: ntrial
    
    [strip.label,strip.trial{i}] = create_bipolar_montage_strip(strip.label,strip.trial{i});
    
    if(isempty(strip.trial{i}))
        strip = [];
        break
    end
    
end

outdata = merge_dataset(gridX,gridY);
outdata = merge_dataset(outdata,strip);
