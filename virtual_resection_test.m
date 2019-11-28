%% load data

addpath(genpath('/home/matteo/Desktop/git_rep/ieeg_respect_bids/'))
addpath('/home/matteo/Desktop/git_rep/virtual_resection/montage/')
addpath('/home/matteo/Desktop/git_rep/fieldtrip/')
addpath('/home/matteo/Desktop/git_rep/epi/matlab/episign/biomarker_pipeline/')
addpath('/home/matteo/Desktop/git_rep/epi/matlab/stiliyan/matLabTools/');

bidsFolder = '/home/matteo/Desktop/tle_e/converted/';
sitFolder  = '/home/matteo/Desktop/tle_e/converted/sub-RESP0465/ses-SITUATION1A/ieeg/';
sitFName   = 'sub-RESP0465_ses-SITUATION1A_task-acute_ieeg';
% read data

[status,msg,cfg,data] = importBidsData(sitFolder);

% use the last minute
cfg         = [];
cfg.trials  = 'all';
cfg.length  = 60; %seconds of new trials
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

%% detrend demean

cfgPre.demean  = 'yes';
cfgPre.detrend = 'yes';
cfgPre.trial   = 'all';
m_data         = ft_preprocessing(cfgPre,m_data);


%% view data
cfg           = [];
cfg.viewmode  = 'vertical';  
cfg.blocksize = 15;
cfg.seldat    = 'current';

ft_databrowser(cfg,m_data)







