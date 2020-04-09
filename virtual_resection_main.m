%  Implementation of the virtual resection with intraoperative intracranial
%  data. 
%  Two types of virtual resection based on functional connectivity are implemented:
%
%   1) Naive virtual resection where resected nodes are removed from the connectivity 
%      matrix (see Kini 2019 Brain) 
%
%   2) Virtual resection with partialization where common source influences 
%      from resected nodes (signals from resected electrodes) are removed before 
%      to compute functional connectivity  
%
%
%  see Demuru et al. doi: 


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


function virtual_resection_main()

% add dependecies
path_settings();

% BIDS input folder 
bidsFolder = '/home/matteo/Desktop/virtual_resection/sel_data/';
% output folder where to save results
outFolder  = '/home/matteo/Desktop/virtual_resection/commenting/';



% subject information table (containing seizure outcome and other info)
info_F = '/home/matteo/Desktop/virtual_resection/info/info.tsv';
info_T = readtable(info_F,'FileType','text','Delimiter','tab','ReadVariableNames',1);
info_T.Properties.RowNames = info_T.subjID;
% subject and situations pre and post to use
subjSit2use_F = '/home/matteo/Desktop/virtual_resection/info/situations2use.tsv';
subjSit2use_T = readtable(subjSit2use_F,'FileType','text','Delimiter','tab','ReadVariableNames',1);
subjSit2use_T.Properties.RowNames = subjSit2use_T.subjID;


% select subjects
subjName = subjSit2use_T.subjID;

%subjName = {'RESP0381','RESP0384','RESP0396','RESP0428','RESP0465','RESP0586','RESP0619','RESP0659'};
%subjName = {'RESP0384'};

  
out            = [];
cfg.bidsFolder = bidsFolder;

cfg.fc_type         = 'h2';  % functional connectivity measure
cfg.montage         = 'avg'; % montage to apply
    
cfg.lastMinute      = 1;     % flag to control the if using the last minute or the whole recordings
cfg.L               = 60;    % last minute
cfg.overlap         = 0;     % no overlapping 

cfg.Tlength         = 5;     % seconds of new trials
cfg.Toverlap        = 0;     % no overlapping between the epochs

cfg.notch           = [50 100 150]; % notch filtering 
cfg.bpfreq          = [5  150];     % frequency interval of interest 

cfg.notchBool       = 'yes';
cfg.bpBool          = 'yes';

cfg.EnvType         = 'filt_hilbert_env'; % compute the envelope using hilbert
cfg.FConEnv         = 'yes';              % flag to control if using envelopes of the signal or raw signal

outFolder = fullfile(outFolder,'syn',cfg.fc_type);
    
if(~exist(outFolder,'dir'))
    mkdir(outFolder)
end
    
for i = 1: numel(subjName)
    
    cfg.subjName        = subjName{i};
    % take the seizure outcome at 1 year
    cfg.seizOut         = info_T{subjName{i},'description_sf_1y'};
    cfg.sitNames        = subjSit2use_T{subjName{i},{'pre','post'}};
    switch cfg.fc_type
        case 'h2'
                out = nonLinear_partialization(cfg);
    end
    % save every situations in a different file 
    for s = 1 : numel(out)
        
        sit_res = out{s};
           
        save(fullfile(outFolder,strcat('sub-',subjName{i},'_','ses-',sit_res{1}.sitName,'_',sit_res{1}.fc_type)),'sit_res');
        
    end
end


% Pre-processing step for the analysis
% notch filtering / frequency band filtering / demean / applying montage  
%
%
% INPUT
%      cfg  - struct with the following fields 
%               bidsFolder: input folder where to find the data in BIDS format (see https://github.com/suforraxi/ieeg_respect_bids) 
%               fc_type:    functional connectivity measure to compute (in our case non-linear correlation 'h2')
%               montage:    reference montage to apply to the data (in our case 'avg')
%               lastMinute: flag to control if perfoming the analysis on the last minute (1) or on the whole data (0)
%               L:          last minute expressed in second (60) (input for fieldtrip routines)
%               overlap:    no overlap flag (0) (input for fieldtrip routines) 
%               Tlength:    epoch length used to perfom the analysis   
%               Toverlap:   no overlap flag (0) 
%               notch:      frequencies to apply notch filter ([50 100 150])
%               bpfreq:     band filtering paramenters  ([5 150])
%               notchBool:  flag to control notch filter (yes/no) 
%               bpBool:     flag to control band pass filetering (yes/no)  
%               EnvType:    define how to compute the envelope ('filt_hilbert_env' using Hilbert transform)
%               FConEnv:    flag to control if compute the functional connectity measure of raw signal or on the envelopes (yes/no)
%               subjName:   selected coded subject name 'RESPXXXX' for whom computing the functional connectivity
%               sitNames:   selected situations per subjects
%               seizOut:    seizure outcome string {i.e. EngelClass_AED_InfoAboutMedication:'1A_AED_stop' }
%
% OUTPUT
%       cfgAnalysis - struct with fields symmarizing the analysis parameters and subject information
%                     sitFName: string with filename of the situation analyze
%                               ('sub-RESPXXXX_ses-SITUATIONXX_task-acute_ieeg')
%                     subjName: string for coded subject name to analyze 
%                               ('RESPXXXX')
%                     sitName:  string with the situation to analyze
%                               ('SITUATIONXX')
%                     bpfreq:   array with [high band pass, low band pass] 
%                               ([5 150])
%                     epoch:    integer representing the epoch length
%                     overlap:  integer used as flag 
%                               to control if use overlapping epochs or not(1/0)
%                     seizOut:  seizure outcome string {i.e. EngelClass_AED_InfoAboutMedication:'1A_AED_stop' }
%                     fc_type:  string representing the functional connectivity measure used'
%                     EnvType:  string representig how the envelope was computed 'filt_hilbert_env'
%                     FConEnv:  string represeing if the functional
%                               connectivity was computed on the envelopes 
%                               or on the raw signals (yes/no)
%                     pre:      integer defining if the situation was pre-resection or not(1/0)
%                     post:     FConEnv: 'yes'
%                     montage:  reference montage applied to the data ('avg' for average /'bipolar' for bipolar in two direction)  
%                     res_ch:   cell array with the name of resected channels
%           
%
%       m_data      - fieldtrip data structure after pre-processing
%
%
function [cfgAnalysis,m_data] = preprocessing_steps(cfg)

% read data
bidsFolder = cfg.bidsFolder; 
subjName   = cfg.subjName;


cfgAnalysis = [];
m_data      = [];


situationName = cfg.situationName{1};
    
if(cfg.pre || cfg.post) % pre or post situation
       
        sitFolder =  fullfile(bidsFolder,strcat('sub-',subjName),strcat('ses-',situationName),'ieeg',filesep);
        sitFName  =  strcat('sub-',subjName,'_','ses-',situationName,'_task-acute_ieeg') ;

        [status,msg,cfgImport,data] = importBidsData(sitFolder);

        
        cfgAnalysis.sitFName = sitFName;
        cfgAnalysis.subjName = subjName;
        cfgAnalysis.sitName  = situationName;
        cfgAnalysis.bpfreq   = cfg.bpfreq;
        cfgAnalysis.epoch    = cfg.Tlength;
        cfgAnalysis.overlap  = cfg.Toverlap;
        cfgAnalysis.seizOut  = cfg.seizOut;
        cfgAnalysis.fc_type  = cfg.fc_type;
        cfgAnalysis.EnvType  = cfg.EnvType;
        cfgAnalysis.pre      = cfg.pre;
        cfgAnalysis.post     = cfg.post;
        cfgAnalysis.FConEnv  = cfg.FConEnv;
        cfgAnalysis.montage  = cfg.montage;
        
        if(status == 0) 
            maxtime = max(cellfun(@length,data.time));
            if(  maxtime > cfg.L*data.fsample )
               
                % select the data to analyze
                
                cfgLastEp         = [];
                cfgLastEp.channel = {'Gr*'};
                
                if(cfg.lastMinute == 1 )
                    
                    cfgRe         = [];
                    cfgRe.trials  = 'all';
                    cfgRe.length  = cfg.L; %seconds of new trials
                    cfgRe.overlap = cfg.overlap;

                    data = ft_redefinetrial(cfgRe,data);
                  
                    % select the grid channels (i.e. Gr)
                   
                    ntrials           = size(data.trial,2);
                    cfgLastEp.trials  = ntrials;
                    
                else % all data available
                    cfgLastEp.trials  = 'all';       
                end
              
                m_data            = ft_selectdata(cfgLastEp,data); 
                
                % apply montage 
                
                switch cfg.montage 
                    
                    case 'avg'    
                            % avg montage
                            cfgM.reref       = 'yes';
                            cfgM.refmethod   = 'avg';
                            cfgM.implicitref = [];
                            cfgM.refchannel  = 'all';
                            m_data           = ft_preprocessing(cfgM,m_data);
                    case 'bipolar'
                     %bipolar two directions
                            m_data = apply_bipolar2D_montage(m_data);
                end
                % remove artefacts

                [ res_channel, artefact_T] = get_metadata(bidsFolder,sitFName);
                
                cfgAnalysis.res_ch         = res_channel;
               
               if (cfg.lastMinute) 
                    
                   [idxChArtefact ,idx_art_trial] = find_artefacts(m_data.sampleinfo,m_data.label,artefact_T);
                    cfgCH.channel                 = m_data.label(~idxChArtefact);
                    m_data                        = ft_preprocessing(cfgCH,m_data);  
               end
                
                % redefine trials
               cfgReTrials.trials  = 'all';
               cfgReTrials.length  = cfg.Tlength; %seconds of new trials
               cfgReTrials.overlap = cfg.Toverlap;

               m_data = ft_redefinetrial(cfgReTrials,m_data);
                 
                %% remove power line (50Hz)
                
               cfgNotch.dftfilter = cfg.notchBool;
               cfgNotch.dftfreq   = cfg.notch;
               cfgNotch.trial    = 'all';
                
               m_data = ft_preprocessing(cfgNotch,m_data);
                
                %% detrend demean
                
               cfgPre.demean   = 'yes';
               cfgPre.detrend  = 'yes';
               cfgPre.trial    = 'all';
               cfgPre.bpfilter = cfg.bpBool;
               cfgPre.bpfreq   = cfg.bpfreq;

               m_data         = ft_preprocessing(cfgPre,m_data);
            end %max time
        end % status
    end % pre post



% compute h2 on the envelope of the signals for pre-resection and post-resection recordings
%
% INPUT 
%       cfg  - struct with the following fields 
%               bidsFolder: input folder where to find the data in BIDS format (see https://github.com/suforraxi/ieeg_respect_bids) 
%               fc_type:    functional connectivity measure to compute (in our case non-linear correlation 'h2')
%               montage:    reference montage to apply to the data (in our case 'avg')
%               lastMinute: flag to control if perfoming the analysis on the last minute (1) or on the whole data (0)
%               L:          last minute expressed in second (60) (input for fieldtrip routines)
%               overlap:    no overlap flag (0) (input for fieldtrip routines) 
%               Tlength:    epoch length used to perfom the analysis   
%               Toverlap:   no overlap flag (0) 
%               notch:      frequencies to apply notch filter ([50 100 150])
%               bpfreq:     band filtering paramenters  ([5 150])
%               notchBool:  flag to control notch filter (yes/no) 
%               bpBool:     flag to control band pass filetering (yes/no)  
%               EnvType:    define how to compute the envelope ('filt_hilbert_env' using Hilbert transform)
%               FConEnv:    flag to control if compute the functional connectity measure of raw signal or on the envelopes (yes/no)
%               subjName:   selected coded subject name 'RESPXXXX' for whom computing the functional connectivity
%               sitNames:   selected situations per subjects
%               seizOut:    seizure outcome string {i.e. EngelClass_AED_InfoAboutMedication:'1A_AED_stop' }
%
%
% OUTPUT
%      out - cell containg the results for a subject. 
%            It contains two cell elements, one for pre-resection situation, one for the post-resection situation.
%            Each element of the cell is a cell summaring the results for
%            every epoch of the situation (see output of compute_fc_h2)
%
function out = nonLinear_partialization(cfg)

out = [];

% input folder and subject name to read the data
 bidsFolder = cfg.bidsFolder; 
 subjName   = cfg.subjName;
  
 sit             = dir(fullfile(bidsFolder,strcat('sub-',subjName),strcat('ses*')));    
 situationName   = cell(numel(cfg.sitNames),1);
 
 % filter selected situation 
 k = 1;
 for i = 1 :numel(sit)
     
     c_situationName = replace(sit(i).name,'ses-','');
    
    if(any(strcmp(c_situationName,cfg.sitNames)))
        situationName{k} = c_situationName;
        k                = k+1;
    end
     
 end

 
 [pre,~,post] = find_pre_int_post(situationName);
 
  
 k = 1;
 for s = 1 : numel(situationName)
     
       cfg.pre           = pre(s);
       cfg.post          = post(s);
       cfg.situationName = situationName(s);

       [cfgAnalysis,m_data] = preprocessing_steps(cfg);
               
       o      = compute_fc_h2(cfgAnalysis,m_data);
       out{k} = o;  
       k      = k + 1;
    
 end

% Compute functional connectivity measure
%
%
% INPUT
%       cfg  - struct with fields summarizing the analysis parameters and subject information
%                     sitFName: string with filename of the situation analyze
%                               ('sub-RESPXXXX_ses-SITUATIONXX_task-acute_ieeg')
%                     subjName: string for coded subject name to analyze 
%                               ('RESPXXXX')
%                     sitName:  string with the situation to analyze
%                               ('SITUATIONXX')
%                     bpfreq:   array with [high band pass, low band pass] 
%                               ([5 150])
%                     epoch:    integer representing the epoch length
%                     overlap:  integer used as flag 
%                               to control if use overlapping epochs or not(1/0)
%                     seizOut:  seizure outcome string {i.e. EngelClass_AED_InfoAboutMedication:'1A_AED_stop' }
%                     fc_type:  string representing the functional connectivity measure used'
%                     EnvType:  string representig how the envelope was computed 'filt_hilbert_env'
%                     FConEnv:  string represeing if the functional
%                               connectivity was computed on the envelopes 
%                               or on the raw signals (yes/no)
%                     pre:      integer defining if the situation was pre-resection or not(1/0)
%                     post:     FConEnv: 'yes'
%                     montage:  reference montage applied to the data ('avg' for average /'bipolar' for bipolar in two direction)  
%                     res_ch:   cell array with the name of resected channels
%           
%
%       m_data - fieldtrip data structure after pre-processing
%
%
%
% OUTPUT 
%
%       o     - cell array, each element is a struct summarizing the results for each epoch. 
%               Each struct contains the following fields
%               
%               fc_type:  string representing the fucntional connectivity measure computed
%               band:     array representing the frequency band interval used to filter the data ([5 150])
%               epoch:    integer representing the epoch length in second 
%               overlap:  0 if there was no overlapping with the epochs 1 with overlapping
%               sitFName: string with the filename of the data input 
%                         according to BIDS format ('sub-RESPXXXX_ses-SITUATIONXX_task-acute_ieeg')
%               subjName: string of the coded subject name from the Utrecht RESPect database ('RESPXXXX')
%               sitName:  string of the situation analyzed 'SITUATIONXX'
%               so:       seizure outcome string {i.e. EngelClass_AED_InfoAboutMedication:'1A_AED_stop' }
%               sitType:  string representing if the situation was pre-resection or post-resection ('Pre'/'Post')
%               C:        2d array representing the functional connectivity matrix 
%               C1:       2d array representing the functional connectivity matrix after the naive virtual resection 
%               C2:       2d array representing the functional connectivity matrix after the virtual resection with partialization 
%

function o = compute_fc_h2(cfg,m_data)



nTrial      = numel(m_data.trial);
res_channel = cfg.res_ch;
fc_type     = cfg.fc_type;
pre         = cfg.pre;
post        = cfg.post;

o      = []; % struct with the result for each epoch of a  situation


for t = 1 : nTrial     

    mEnv            = get_Envelope(cfg,m_data.trial{t});
    C               = fc(mEnv,fc_type);
    
 
    o{t}.fc_type  = fc_type;
    o{t}.band     = cfg.bpfreq; 
    o{t}.epoch    = cfg.epoch;
    o{t}.overlap  = cfg.overlap;
  
    o{t}.C        = C;
  
    o{t}.sitFName = cfg.sitFName;
    o{t}.subjName = cfg.subjName;
    o{t}.sitName  = cfg.sitName;
    o{t}.so       = cfg.seizOut;
    o{t}.montage  = cfg.montage;

end
for t = 1 : nTrial % virtual resection
    if(pre)
        o{t}.sitType = 'Pre';

         % virtual resection
        if(numel(res_channel) > 0) % there are resected channels

            % virtual resection naive
            idx2rm = zeros(numel(m_data.label),1);
            % find indexes to be resected
            for r = 1 :numel(res_channel)
                 mono_ch = res_channel{r};
                 idx2rm  = idx2rm | (~cellfun(@isempty,regexp(m_data.label,mono_ch)));
            end
            C           = o{t}.C;
            C1          = C(~idx2rm,~idx2rm);
            o{t}.C1     = C1;
        
            % virtual resection with partialization
            m_data   = removeNoNLinear_ch(m_data,idx2rm); 
            mEnv     = get_Envelope(cfg,m_data.trial{t});
            C2       = fc(mEnv,fc_type);
            o{t}.C2  = C2;
            
        else % no need to remove channels (it should not happen given the data)
            o{t}.C1     = C;
            o{t}.C2     = C;
        end
    elseif(post)
        o{t}.sitType = 'Post';
        o{t}.C1      = NaN;
        o{t}.C2      = NaN;
       
    end
end


