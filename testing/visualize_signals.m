function visualize_signals(sitFolder)

% /home/matteo/Desktop/tle_e/converted/sub-RESP0634/ses-SITUATION3A/ieeg
bidsFolder = regexp(sitFolder,'sub-\w*','split');
bidsFolder = bidsFolder{1};
aux        = regexp(sitFolder,'\w*(?<subName>sub-\w*).(?<sesName>ses-\w*).ieeg','names');     
sitFName   = strcat(aux.subName,'_',aux.sesName,'_','task-acute_ieeg'); 



L = 5; % last minute
[status,msg,cfg,data] = importBidsData(sitFolder);

maxtime = max(cellfun(@length,data.time));

%if(status == 0 && maxtime > L*data.fsample )
if(status == 0  )

    
            cfg         = [];
            cfg.trials  = 'all';
            cfg.length  = L; %seconds of new trials
            cfg.overlap = 0;

            data = ft_redefinetrial(cfg,data);

            %cfgLastEp        = [];
            %ntrials          = size(data.trial,2);
            %cfgLastEp.trials = ntrials;
            %data             = ft_selectdata(cfgLastEp,data); 

            % apply montage 

            %m_data = apply_bipolar2D_montage(data);

            % avg montage
            cfgM.reref       = 'yes';
            cfgM.refmethod   = 'avg';
            cfgM.implicitref = [];
            cfgM.refchannel  = 'all';
            m_data           = ft_preprocessing(cfgM,data);

            % remove artefacts

            [ res_channel, artefact_T]     = get_metadata(bidsFolder,sitFName);
            [idxChArtefact ,idx_art_trial] = find_artefacts(m_data.sampleinfo,m_data.label,artefact_T);
            cfgCH.channel                  = m_data.label(~idxChArtefact);
            m_data                         = ft_preprocessing(cfgCH,m_data);  

            %% detrend demean

            cfgPre.demean  = 'yes';
            cfgPre.detrend = 'yes';
            cfgPre.trial   = 'all';
            m_data         = ft_preprocessing(cfgPre,m_data);


    cfg           = [];
    cfg.viewmode  = 'vertical';  
    cfg.blocksize = 15;
    cfg.seldat    = 'current';

    ft_databrowser(cfg,m_data)

else
    'not enough time'
end