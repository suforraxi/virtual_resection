% load data

function coh_before_after_main()


addpath('/home/matteo/Desktop/git_rep/epi/matlab/episign/biomarker_pipeline/')
addpath('/home/matteo/Desktop/git_rep/epi/matlab/stiliyan/matLabTools/');

bidsFolder = '/home/matteo/Desktop/virtual_resection/sel_data/';
outFolder  = '/home/matteo/Desktop/virtual_resection/';



info_F = '/home/matteo/Desktop/virtual_resection/info/info.tsv';
info_T = readtable(info_F,'FileType','text','Delimiter','tab','ReadVariableNames',1);
info_T.Properties.RowNames = info_T.subjID;


subjName = {'RESP0381','RESP0384','RESP0396','RESP0428','RESP0465','RESP0586','RESP0619','RESP0659'};
%subjName = {'RESP0381','RESP0451'};
  
out            = [];
cfg.bidsFolder = bidsFolder;

cfg.fc_type         = 'h2';
cfg.montage         = 'avg';
    
cfg.lastMinute      = 1;
cfg.L               = 60;
cfg.overlap         = 0;

cfg.Tlength         = 5; %seconds of new trials
cfg.Toverlap        = 0;

cfg.notch           = [50 100 150];
cfg.bpfreq          = [5  150];

cfg.notchBool       = 'yes';
cfg.bpBool          = 'yes';

cfg.EnvType         = 'filt_hilbert_env';
cfg.FConEnv         = 'yes';


outFolder = fullfile(outFolder,'syn',cfg.fc_type);
    
if(~exist(outFolder,'dir'))
    mkdir(outFolder)
end
    
for i = 1: numel(subjName)
    
    cfg.subjName        = subjName{i};
    cfg.seizOut         = info_T{subjName{i},'description_sf_1y'};
 
    switch cfg.fc_type
        case {'corr','coh','pli'}
                out = linear_partialization(cfg);
        case 'h2'
                out = nonLinear_partialization(cfg);
    end
  
    for s = 1 : numel(out)
        
        sit_res = out{s};
           
        save(fullfile(outFolder,strcat('sub-',subjName{i},'_','ses-',sit_res{1}.sitName,'_',sit_res{1}.fc_type)),'sit_res');
        
    end
end


function out = linear_partialization(cfg)

out = [];
%L   = cfg.L; % time in second for the epoch

% read data
 bidsFolder = cfg.bidsFolder; 
 subjName   = cfg.subjName;
 seizOut    = cfg.seizOut;

 sit             = dir(fullfile(bidsFolder,strcat('sub-',subjName),strcat('ses*')));    
 situationName   = cell(numel(sit),1);
 
 for i = 1 :numel(sit)
     situationName{i} = replace(sit(i).name,'ses-','');
 end

[pre,inter,post] = find_pre_int_post(situationName);

k = 1;
for s = 1 : numel(situationName)
     
       cfg.pre           = pre(s);
       cfg.post          = post(s);
       cfg.situationName = situationName(s);
  
       [cfgAnalysis,m_data] = preprocessing_steps(cfg);
                
       o = compute_fc(cfgAnalysis,m_data);
                
       out{k} = o;  
       k      = k + 1;
    
end



function [cfgAnalysis,m_data] = preprocessing_steps(cfg)

% read data
bidsFolder = cfg.bidsFolder; 
subjName   = cfg.subjName;
seizOut    = cfg.seizOut;

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





function out = nonLinear_partialization(cfg)

out = [];
%L   = cfg.L; % time in second for the epoch

% read data
 bidsFolder = cfg.bidsFolder; 
 subjName   = cfg.subjName;
 seizOut    = cfg.seizOut;
% 
 sit             = dir(fullfile(bidsFolder,strcat('sub-',subjName),strcat('ses*')));    
 situationName   = cell(numel(sit),1);
 
 for i = 1 :numel(sit)
     situationName{i} = replace(sit(i).name,'ses-','');
 end
% 
 [pre,inter,post] = find_pre_int_post(situationName);
% 
  
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


function o = compute_fc(cfg,m_data)
 
nTrial      = numel(m_data.trial);
res_channel = cfg.res_ch;
fc_type     = cfg.fc_type;
pre         = cfg.pre;
post        = cfg.post;
o           = []; % struct with the result for situation
              

for t = 1 : nTrial
% compute functional connectivity measure on envelopes
% or on the raw signals
    %L = [];
    C = [];
    if(strcmp(cfg.FConEnv,'yes'))

        mEnv         = get_Envelope(cfg,m_data.trial{t});
        C            = fc(mEnv,fc_type);
    else% kini 2019 no envelope 
        C            = fc(m_data.trial{t},fc_type);
    end
    L            = get_laplacian(C);
    [syn, V, E]  = get_synchronizability(L);
    
    o{t}.FConEnv  = cfg.FConEnv;
    o{t}.fc_type  = cfg.fc_type;
    o{t}.band     = cfg.bpfreq; 
    o{t}.epoch    = cfg.epoch;
    o{t}.overlap  = cfg.overlap;
    o{t}.syn      = syn;
    o{t}.C        = C;
    o{t}.L        = L;
    o{t}.V        = V;
    o{t}.E        = E;
    o{t}.sitFName = cfg.sitFName;
    o{t}.subjName = cfg.subjName; 
    o{t}.sitName  = cfg.sitName;
    o{t}.so       = cfg.seizOut;

end
for t = 1 : nTrial % virtual resection
    if(pre)
        o{t}.sitType = 'Pre';

         % virtual resection
        if(numel(res_channel) > 0) % there are resected channels

            % virtual resection without orthogonalization
            idx2rm = zeros(numel(m_data.label),1);
            for r = 1 :numel(res_channel)
                 mono_ch = res_channel{r};
                 idx2rm  = idx2rm | (~cellfun(@isempty,regexp(m_data.label,mono_ch)));

            end
            C           = o{t}.C;
            C1          = C(~idx2rm,~idx2rm);
            L1          = get_laplacian(C1);
            [syn, V, E] = get_synchronizability(L1);
            o{t}.V1syn  = syn;
            o{t}.C1     = C1;
            o{t}.L1     = L1;
            o{t}.V1     = V;
            o{t}.E1     = E;   

            for r = 1 :numel(res_channel)
                 % virtual resection with orthogonalization
                 mono_ch = res_channel{r};
                 idx2rm  = find(~cellfun(@isempty,regexp(m_data.label,mono_ch)));
                while(~isempty(idx2rm))
                      m_data = remove_ch(m_data,idx2rm(1));
                      idx2rm = find(~cellfun(@isempty,regexp(m_data.label,mono_ch)));
                end
                
                C2 = [];
                
                if(strcmp(cfg.FConEnv,'yes'))
                    mEnv = get_Envelope(cfg,m_data.trial{t});
                    C2   = fc(mEnv,fc_type);                
                else
                    C2   = fc(m_data.trial{t},fc_type);    
                end
                
                L2              = get_laplacian(C2);
                [syn, V, E]     = get_synchronizability(L2);
                o{t}.V2syn      = syn;
                o{t}.C2         = C2;
                o{t}.L2         = L2;
                o{t}.V2         = V;
                o{t}.E2         = E;

            end

        end 
    elseif(cfg.post)
        o{t}.sitType = 'Post';

        o{t}.V1syn     = NaN;
        o{t}.C1        = NaN;
        o{t}.L1        = NaN;
        o{t}.V1        = NaN;
        o{t}.E1        = NaN;

        o{t}.V2syn      = NaN;
        o{t}.C2         = NaN;
        o{t}.L2         = NaN;
        o{t}.V2         = NaN;
        o{t}.E2         = NaN;
    end
end



% Compute three things:
% 
% 1) Functional connectivity pre-resection and post-resection.
%
% 2) Simulate the resection from pre-resection recordings 
%    (virtual resection naive) as in the paper of Kini 2019 
%    (naive removing the channels without removing 'the effect')
%
% 3) Simulate the resection from pre-resection recordings
%    (virtual resection with signal partialization) 
%
%

function o = compute_fc_h2(cfg,m_data)



nTrial      = numel(m_data.trial);
res_channel = cfg.res_ch;
fc_type     = cfg.fc_type;
pre         = cfg.pre;
post        = cfg.post;

o      = []; % struct with the result for situation


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


