% load data

function coh_before_after_main()


addpath('/home/matteo/Desktop/git_rep/epi/matlab/episign/biomarker_pipeline/')
addpath('/home/matteo/Desktop/git_rep/epi/matlab/stiliyan/matLabTools/');

bidsFolder = '/home/matteo/Desktop/virtual_resection/sel_data/';
outFolder  = '/home/matteo/Desktop/virtual_resection/';



info_F = '/home/matteo/Desktop/virtual_resection/info/info.tsv';
info_T = readtable(info_F,'FileType','text','Delimiter','tab','ReadVariableNames',1);
info_T.Properties.RowNames = info_T.subjID;


subjName = {'RESP0381','RESP0384','RESP0428','RESP0451','RESP0465','RESP0586','RESP0659'};
%subjName = {'RESP0381'};
        
cfg.bidsFolder =  bidsFolder;
out            = [];
rTbl           = [];

for i = 1: numel(subjName)
    
    cfg.subjName        = subjName{i};
    cfg.seizOut         = info_T{subjName{i},'description_sf_1y'};
    cfg.fc_type         = 'corr';
    cfg.L               = 60;
    cfg.overlap         = 0;
    
    cfg.Tlength         = 5; %seconds of new trials
    cfg.Toverlap        = 0;
    
    cfg.notch           = [50 100 150];
    cfg.bpfreq          = [5   150];
    %cfg.bpfreq          = [0   0];
   
    cfg.notchBool       = 'yes';
    cfg.bpBool          = 'yes';
    
    
    cfg.type            = 'filt_hilbert_env';
    %cfg.type            = 'huang_env';
    cfg.nHComp          = 3;
    cfg.selected_comp   = 2;
    
    
    %cfg.NbestTrial      = 4;
    %cfg.ScoreFreqRange  = 5:100;
    %cfg.ScoreSmoothF    = 3;
    
    
    out                 = coh_before_after(cfg);
    
    
    for s = 1 : numel(out)
        
        sit_res = out{s};
        
        %save(fullfile(outFolder,'coh',strcat('sub-',subjName{i},'_','ses-',sit_res{1}.sitName,'_',sit_res{1}.fc_type)),'sit_res');
        save(fullfile(outFolder,'coh',strcat('sub-',subjName{i},'_','ses-',sit_res{1}.sitName)),'sit_res');
        
        %rTbl = [rTbl; create_table(sit_res) ];
        
    end
end

%save(fullfile(outFolder,'tables','summary_tbl_coh'),'rTbl');

function out = coh_before_after(cfg)

out = [];
L   = cfg.L; % time in second for the epoch

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

fc_type = cfg.fc_type;
k = 1;
for s = 1 : numel(situationName)
    
    if(pre(s) || post(s)) % pre or post situation
       
        sitFolder =  fullfile(bidsFolder,strcat('sub-',subjName),strcat('ses-',situationName{s}),'ieeg',filesep);
        sitFName  =  strcat('sub-',subjName,'_','ses-',situationName{s},'_task-acute_ieeg') ;

        [status,msg,cfgImport,data] = importBidsData(sitFolder);

        
        if(status == 0) 
            maxtime = max(cellfun(@length,data.time));
            if(  maxtime > L*data.fsample )
                % use the last minute
                cfgRe         = [];
                cfgRe.trials  = 'all';
                cfgRe.length  = L; %seconds of new trials
                cfgRe.overlap = cfg.overlap;

                data = ft_redefinetrial(cfgRe,data);

                cfgLastEp         = [];
                ntrials           = size(data.trial,2);
                cfgLastEp.trials  = ntrials;
                cfgLastEp.channel = {'Gr*'};
                m_data            = ft_selectdata(cfgLastEp,data); 
                
                
                
                % apply montage 
                % avg
                cfgM.reref       = 'yes';
                cfgM.refmethod   = 'avg';
                cfgM.implicitref = [];
                cfgM.refchannel  = 'all';
                m_data           = ft_preprocessing(cfgM,m_data);
                
                % bipolar two directions
                %m_data = apply_bipolar2D_montage(m_data);

                % remove artefacts

                [ res_channel, artefact_T]     = get_metadata(bidsFolder,sitFName);
                [idxChArtefact ,idx_art_trial] = find_artefacts(m_data.sampleinfo,m_data.label,artefact_T);
                cfgCH.channel                  = m_data.label(~idxChArtefact);
                m_data                         = ft_preprocessing(cfgCH,m_data);  
                
                %% Select 'homogeneous trials'
               % cfgScore.freqRange      = cfg.ScoreFreqRange;
               % cfgScore.fs             = int32(m_data.fsample);
               % cfgScore.windowL        = cfg.Tlength;
               % cfgScore.smoothFactor   = cfg.ScoreSmoothF;
                
                %[idx_best_trial,~,scores] = scorEpochs(cfgScore,m_data.trial{1});
                
                %% redefine trials
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
                

                %cfgSelT.trials  = idx_best_trial(1:cfg.NbestTrial);
                %m_data          = ft_selectdata(cfgSelT,m_data);
                
                
                % compute envelope and coherence
                nTrial = numel(m_data.trial);
                o      = []; % struct with the result for situation
                for t = 1 : nTrial
                %if(nTrial == 1 )      

                    mEnv            = get_Envelope(cfg,m_data.trial{t});
                    C               = fc(mEnv,fc_type);
                    m_data.trial{t} = mEnv;
                    [coh, V, E]     = get_coherence(C);
                    
                    
                    o{t}.fc_type  = fc_type;
                    o{t}.band     = cfgPre.bpfreq; 
                    o{t}.epoch    = cfgReTrials.length;
                    o{t}.overlap  = cfgReTrials.overlap;
                    o{t}.coh      = coh;
                    o{t}.C        = C;
                    o{t}.V        = V;
                    o{t}.E        = E;
                    o{t}.sitFName = sitFName;
                    o{t}.subjName = subjName; 
                    o{t}.sitName  = situationName{s};
                    o{t}.so       = seizOut;
                    
                    
                    
                    
                end
                for t = 1 : nTrial % virtual resection
                    if(pre(s))
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
                            C1           = C(~idx2rm,~idx2rm);
                            [coh, V, E] = get_coherence(C1);
                            o{t}.V1coh  = coh;
                            o{t}.C1     = C1;
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
                                C2              = fc(m_data.trial{t},fc_type);
                                [coh, V, E]     = get_coherence(C2);
                                o{t}.V2coh      = coh;
                                o{t}.C2         = C2;
                                o{t}.V2         = V;
                                o{t}.E2         = E;
                                
                            end
                                            
                        end 
                    elseif(post(s))
                        o{t}.sitType = 'Post';
                       
                        o{t}.V1coh     = NaN;
                        o{t}.C1        = NaN;
                        o{t}.V1        = NaN;
                        o{t}.E1        = NaN;
                        
                        o{t}.V2coh      = NaN;
                        o{t}.C2         = NaN;
                        o{t}.V2         = NaN;
                        o{t}.E2         = NaN;
                    end
                end
                
                
                
                out{k} = o;  
                k      = k + 1;
        end % max time
    end%status
end % pre or post
    
    
end


% create result table
%
% INPUT
%  res   - struct with the following fields
%                 res.coh      - coherence computed as the ratio between the maximum of weighted eigenvalue and the sum of the weighted eigenvalues    
%                 res.C        - functional connectivity matrix (channel X channel)  
%                 res.V        - eigenvector matrix computed on functional connectivity matrix (channel X #eigenvalues)   
%                 res.E        - eigenvalues matrix computed on the functional connectivity matrix  
%                 res.sitFName - situation file name (sub-RESPXXXX_ses-SITUATIONXX_task-acute_ieeg)
%                 res.subjName - subject name (RESPXXXX)
%                 res.sitName  - situation name (SITUATIONXX)
%                 res.sitType  - string related to the situation type (Pre for pre-resection, Post for post-resection) 
%                 res.so      - string for seizure outcome (i.e. 1A_AED_stop / 1A_AED_low / 1A_AED_eq )  
% OUTPUT
%    
%  r_table - table with the following variables
%                  subjName  - subject name (RESPXXXX)
%                  sitName   - situation name (SITUATIONXX)
%                  sitType   - string related to the situation type (Pre for pre-resection, Post for post-resection) 
%                  coh       - coherence computed as the ratio between the maximum of weighted eigenvalue and the sum of the weighted eigenvalues  
%                  seizOout  - seizure outcome for the patient

function r_tbl  = create_table(res)


         r_tbl  = table( {res.subjName}, {res.sitName}, {res.sitType}, [res.coh], [res.Vcoh] ,[res.so],...
                        'VariableNames',{'subjName','sitName','sitType','biomarker','v_biomarker','seizOut'} ...
                       );

