% load data

function coh_before_after_main()


addpath('/home/matteo/Desktop/git_rep/epi/matlab/episign/biomarker_pipeline/')
addpath('/home/matteo/Desktop/git_rep/epi/matlab/stiliyan/matLabTools/');

bidsFolder = '/home/matteo/Desktop/virtual_resection/sel_data/';
outFolder  = '/home/matteo/Desktop/virtual_resection/';



info_F = '/home/matteo/Desktop/virtual_resection/info/info.tsv';
info_T = readtable(info_F,'FileType','text','Delimiter','tab','ReadVariableNames',1);
info_T.Properties.RowNames = info_T.subjID;

%subjName = {'RESP0320','RESP0384','RESP0451','RESP0619','RESP0634','RESP0168','RESP0288','RESP0591','RESP0659'};

% subjName = {'RESP0067','RESP0118','RESP0124','RESP0149','RESP0168','RESP0231','RESP0288',...
%             'RESP0301','RESP0362','RESP0381','RESP0384','RESP0396','RESP0409','RESP0424',...
%             'RESP0428','RESP0451','RESP0465','RESP0586','RESP0591','RESP0619','RESP0634',...
%             'RESP0659'
%             };

subjName = {'RESP0381', 'RESP0659'};

        
cfg.bidsFolder =  bidsFolder;
out            = [];
rTbl           = [];

for i = 1: numel(subjName)
    
    cfg.subjName = subjName{i};
    cfg.seizOut  = info_T{subjName{i},'description_sf_1y'};
    out          = coh_before_after(cfg);
    
    for s = 1 : numel(out)
        
        sit_res = out{s};
        
        save(fullfile(outFolder,'coh',strcat('sub-',subjName{i},'_','ses-',sit_res.sitName)),'sit_res');
        
        rTbl = [rTbl; create_table(sit_res) ];
        
    end
end

save(fullfile(outFolder,'tables','summary_tbl_coh'),'rTbl');

function out = coh_before_after(cfg)

out = [];
L   = 60; % time in second for the epoch

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
    
    if(pre(s) || post(s)) % pre or post situation
       
        sitFolder =  fullfile(bidsFolder,strcat('sub-',subjName),strcat('ses-',situationName{s}),'ieeg',filesep);
        sitFName  =  strcat('sub-',subjName,'_','ses-',situationName{s},'_task-acute_ieeg') ;

        [status,msg,cfg,data] = importBidsData(sitFolder);

        
        if(status == 0) 
            maxtime = max(cellfun(@length,data.time));
            if(  maxtime > L*data.fsample )
                % use the last minute
                cfg         = [];
                cfg.trials  = 'all';
                cfg.length  = L; %seconds of new trials
                cfg.overlap = 0;

                data = ft_redefinetrial(cfg,data);

                cfgLastEp        = [];
                ntrials          = size(data.trial,2);
                cfgLastEp.trials = ntrials;
                data             = ft_selectdata(cfgLastEp,data); 

                % apply montage 

                m_data = apply_bipolar2D_montage(data);

                % remove artefacts

                [ res_channel, artefact_T]     = get_metadata(bidsFolder,sitFName);
                [idxChArtefact ,idx_art_trial] = find_artefacts(m_data.sampleinfo,m_data.label,artefact_T);
                cfgCH.channel                  = m_data.label(~idxChArtefact);
                m_data                         = ft_preprocessing(cfgCH,m_data);  
                
                %% redefine trials
%                 cfgReTrials.trials  = 'all';
%                 cfgReTrials.length  = 5; %seconds of new trials
%                 cfgReTrials.overlap = 0;
% 
%                 m_data = ft_redefinetrial(cfgReTrials,m_data);
%                 
                
                %% detrend demean
                
                cfgPre.demean   = 'yes';
                cfgPre.detrend  = 'yes';
                cfgPre.trial    = 'all';
                cfgPre.bpfilter = 'yes';
                cfgPre.bpfreq   = [20 120];

                m_data         = ft_preprocessing(cfgPre,m_data);


                % compute envelope and coherence
                nTrial = numel(m_data.trial);
                o      = []; % struct with the result for situation

                if(nTrial == 1 )      

                    mEnv            = get_Envelope(m_data.trial{1});
                    C               = fc(mEnv);
                    m_data.trial{1} = mEnv;
                    [coh, V, E]     = get_coherence(C);

                    o.coh      = coh;
                    o.C        = C;
                    o.V        = V;
                    o.E        = E;
                    o.sitFName = sitFName;
                    o.subjName = subjName; 
                    o.sitName  = situationName{s};
                    o.so       = seizOut;
                    if(pre(s))
                        o.sitType = 'Pre';
                        % virtual resection
                        if(numel(res_channel) > 0) % there are resected channels
                            
                            for r = 1 :numel(res_channel)
                                mono_ch = res_channel{r};
                                 idx2rm  = find(~cellfun(@isempty,regexp(m_data.label,mono_ch)));
                                while(~isempty(idx2rm))
                                     m_data = remove_ch(m_data,idx2rm(1));
                                     idx2rm = find(~cellfun(@isempty,regexp(m_data.label,mono_ch)));
                                end
                                C           = fc(m_data.trial{1});
                                [coh, V, E] = get_coherence(C);
                                o.Vcoh      = coh;
                                o.VC        = C;
                                o.VV        = V;
                                o.VE        = E;
                                
                            end
                                            
                        end 
                    elseif(post(s))
                        o.sitType = 'Post';
                       
                        o.Vcoh      = NaN;
                        o.VC        = NaN;
                        o.VV        = NaN;
                        o.VE        = NaN;
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

