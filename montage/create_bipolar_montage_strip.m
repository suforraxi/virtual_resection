%create bipolar montage for strip and bids

% ch_label - channel names 
% data     - (ch X samples) data matrix 

% strip2use - strip labels available for which it was possible to compute
%             the bipolar montage
% outdata   - bipolar trasformation for the strip (i.e Str1-Str2)
function [strips2use,outdata]=create_bipolar_montage_strip(ch_label,data)


addhyphen = @(x) [x '-'];

strip2use = {};
outdata   = [];

% look for strip names
strip_names  = {'Str','Rst','Riv','St','StrB','Sst','S1-'};

lookup4strip = regexpi(ch_label,'(Str|Rst|Riv|St|StrB|Sst|S1-)\d+');
strip_idx    = ~cell2mat(cellfun(@isempty,lookup4strip,'UniformOutput',false));
if(any(strip_idx))

    % extract only the strips
    strip_label  = ch_label(strip_idx);
    data         = data(strip_idx,:);

    % order the strips
    [~,reidx]    = sort(strip_label);
    strip_label  = strip_label(reidx);
    data         = data(reidx,:);

    % extract the prefix to construct the montage
    strip_prefix  = regexp(strip_label,'(?<stripname>\w+\D?[0]?)\d+','names');
    unique_prefix = {};

    for i = 1:numel(strip_prefix)
        unique_prefix{i} = strip_prefix{i}.stripname;
    end

    unique_prefix = unique(unique_prefix);
    
    %if(numel(unique_prefix)>1)
    %    error('more than one strip')
    %end
    
    for k = 1 : numel(unique_prefix)
    

        montage = { [1 2] ;...
                    [2 3] ;...
                    [3 4] ;...
                    [4 5] ;...
                    [5 6] ;...
                    [6 7] ;...
                    [7 8] ;...
                  };

        montage_M = zeros(size(montage,1),size(data,1));

            for i = 1:size(montage,1)

                c_mont      = montage{i};
                ch2check    = cell(size(c_mont));
                chAvailable = zeros(numel(strip_label),1);

                for j = 1:numel(c_mont)

                    ch2check{j} = sprintf( '%s%i' , unique_prefix{k}, c_mont(j) );

                    chAvailable = chAvailable | strcmp( strip_label ,ch2check{j});

                end

                curr_ch_idx = strcmp( strip_label , ch2check{1});

                if(sum(chAvailable) == numel(c_mont)) %all required channels are present

                    %curr_ch_idx                        = strcmp( strip_label , ch2check{1}); 
                    chAvailable(curr_ch_idx) = 0;
                    montage_M(i,chAvailable) = -1;
                    montage_M(i,curr_ch_idx) = 1; 

                    %label2add = cellfun(addhyphen,strip_label,'UniformOutput',false);
                    strip2use{i,1} = strcat(ch2check{1},'-',ch2check{end});
                else
                    montage_M(i,:) = 0;
                    strip2use{i,1} = strcat(ch2check{1},'N-N',ch2check{end});
                end

            end


            montage2use = montage_M;
            
          
            outdata{k}    = montage2use * data;
            strips2use{k} = strip2use;
    end
else
    strips2use = {};
    outdata   = [];
end


for k = 1 : size(strips2use,2)
    c_strip  = strips2use{k};
    label    = regexpi(c_strip,'(?<startStr>(Str|Rst|Riv|St|StrB|Sst|S1-))(?<first>\d*)(?<N1>N?)-(?<N2>N?)(?<endStr>(Str|Rst|Riv|St|StrB|Sst|S1-))(?<second>\d*)','names');
    twoDigit = @(x) {sprintf('%02.0f',str2num(x.first)); sprintf('%02.0f',str2num(x.second))};

    padded   = cellfun(twoDigit,label,'UniformOutput',false);
    new_label = cell(size(c_strip));
    for i = 1 : numel(label)

        new_label{i} = strcat(label{i}.startStr,padded{i}{1},label{i}.N1,'-',label{i}.N2,label{i}.endStr,padded{i}{end}); 
    end
    strips2use{k} = new_label;
end