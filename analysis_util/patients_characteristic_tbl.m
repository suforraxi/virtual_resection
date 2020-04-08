% creating a summary table for patient info

% name / age / gender / # electrodes / # resected electrodes / primary
% pathology / 1y seizure outcome


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


function patients_characteristic_tbl()


age_gender_F = '/home/matteo/Desktop/openclinica/TAB_age_gender_second_attempt_2019-10-28-095433974.tsv';
info_F       = '/home/matteo/Desktop/virtual_resection/info/info.tsv';

outFile      = '/home/matteo/Desktop/virtual_resection/info/characteristic_table.txt';
% subject of interest
subj_OI      = {'RESP0381','RESP0384','RESP0396','RESP0428','RESP0465','RESP0586','RESP0619','RESP0659'};

age_gender_T = readtable(age_gender_F,'FileType','text','Delimiter','tab','ReadVariableNames',1);
info_T       = readtable(info_F,'FileType','text','Delimiter','tab','ReadVariableNames',1);

primary_path = {'High Grade Tumor (WHO III + IV)',...
                    'Low Grade Tumor (WHO I + II)',   ...
                    'MTS',                            ...
                    'FCD',                            ...
                    'no abnormalities',               ...
                    'cavernoma',                      ...
                    'gliosis/scar',                   ...
                    'AVM',                            ...
                    'malformation cort. development', ...
                    'TuberoSclerosis'
                };
            
            


% select variables of interest
info_T = info_T(:,{'subjID','description_sf_1y','primary_path_class'});
info_T.Properties.RowNames = info_T.subjID;

% filter for subject of interest
info_T = info_T(subj_OI,:);


age_gender_T = age_gender_T(:,{'StudySubjectID','Sex','DateOfBirth'});
age_gender_T.Properties.VariableNames = {'subjID','Gender','Age'};

current_year  = datetime(date).Year;
dateOFbirth   = age_gender_T.Age;
dateOFbirth   = dateOFbirth.Year;
newAge        = repmat(current_year,numel(dateOFbirth),1) - dateOFbirth;

age_gender_T.Age =  newAge;



characteristics_T = innerjoin(info_T,age_gender_T,'Keys','subjID');

characteristics_T;


infoRec_T = get_info_recordings();

characteristics_T = innerjoin(characteristics_T,infoRec_T,'Keys','subjID');


primary_path_str = cell(numel(subj_OI),1);

for i = 1 : numel(primary_path)
    
    idx = characteristics_T.primary_path_class == i;
    if(any(idx))
        primary_path_str(idx) = primary_path(i);
    end
end

primary_path_str;

final_T = [characteristics_T(:,[1 2 4:10]) cell2table(primary_path_str,'VariableNames',{'primary_pathology'}) ];

writetable(final_T,outFile,'FileType','text', 'Delimiter','tab')

% extract information pre and post resection about:
% number of epochs / number of channels / number of resected channels

function infoRec_T =get_info_recordings()

% folder with the 
inFolder = '/home/matteo/Desktop/virtual_resection/coh/h2_all_epochs/';
subj_OI  = {'RESP0381','RESP0384','RESP0396','RESP0428','RESP0465','RESP0586','RESP0619','RESP0659'};
varNames = {'n_ch_pre','n_RESch_pre','n_Ep_pre','n_ch_post','n_Ep_post'};

m_info = zeros(numel(subj_OI),numel(varNames));

for s = 1 : numel(subj_OI)
    
    pre_post_F = dir(fullfile(inFolder,strcat('*',subj_OI{s},'*')));
    
    for i = 1 : numel(pre_post_F)
        
        % load results
        load(fullfile(pre_post_F(i).folder,pre_post_F(i).name));
        
        nCh_pre     = 0;
        nRESCh_pre  = 0;
        nCh_post    = 0;
       
        nEp_pre     = 0;
        nEp_post    = 0;
        
        switch sit_res{1}.sitType
            case 'Pre'
                  nCh_pre    = size(sit_res{1}.C,1);
                  nRESCh_pre = nCh_pre - size(sit_res{1}.C2,1);
                  nEp_pre    = numel(sit_res);
                  
                  m_info(s,1:3)   = [nCh_pre nRESCh_pre nEp_pre];
            
            case 'Post'
                  nCh_post   = size(sit_res{1}.C,1);  
                  nEp_post   = numel(sit_res);
                  
                  m_info(s,4:end) = [nCh_post nEp_post];
        end
        
    
    end
    
    
end

infoRec_T = array2table(m_info,'VariableNames',varNames);
infoRec_T = [cell2table(subj_OI','VariableNames',{'subjID'}) infoRec_T ];






