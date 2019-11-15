%% load data

addpath(genpath('/home/matteo/Desktop/git_rep/ieeg_respect_bids/'))
addpath('/home/matteo/Desktop/git_rep/virtual_resection/montage/')
addpath('/home/matteo/Desktop/git_rep/fieldtrip/')
addpath('/home/matteo/Desktop/git_rep/epi/matlab/episign/biomarker_pipeline/')
addpath('/home/matteo/Desktop/git_rep/epi/matlab/stiliyan/matLabTools/');

bidsFolder = '/home/matteo/Desktop/tle_e/converted/';
sitFolder  = '/home/matteo/Desktop/tle_e/converted/sub-RESP0465/ses-SITUATION1A/ieeg/';
sitFname   = 'sub-RESP0465_ses-SITUATION1A_task-acute_ieeg';
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

% montage 

% grid
mX_data    = apply_montage2data(data,@create_bipolar_montage_5by4X);
mY_data    = apply_montage2data(data,@create_bipolar_montage_5by4Y);
grid_data  = merge_dataset(mX_data,mY_data);

% strip

strip_data = apply_montage2data_strip(data,@create_bipolar_montage_strip);
m_data     = merge_dataset(grid_data,strip_data);

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

%cfg=[];
%cfg.refchannel={'all'}; 
%m_data = ft_preprocessing(cfg,m_data);

%% ICA or PCA
cfg              = [];
cfg.method       = 'fastica';
cfg.channel      = 'all';
cfg.trials       = 'all';
cfg.numcomponent = 'all';
cfg.demean       = 'yes';
cfg.updatesens   = 'yes';
cfg.feedback     = 'text';


comp = ft_componentanalysis(cfg,m_data);

%% plot components


%% browse ICA

 cfg           = [];
 cfg.layout    = 'grid_5x4_lay_square.mat';
 cfg.zlim      = 'maxabs';
 cfg.highlight = 'labels';
 ft_icabrowser(cfg,comp)
 

%% plot compenent and time-course
cfg            = [];
cfg.layout     = 'grid_5x4_lay_square.mat';
cfg.viewmode   = 'component';
cfg.continuous = 'yes';
cfg.blocksize  = 60;
cfg.channels   = 'all';
cfg.zlim       = 'maxabs';
ft_databrowser(cfg,comp);

%% topoplot of components
figure
for i = 1 : numel(comp.label)
    
    subplot(5,4,i)
    
    cfg           = [];
    cfg.layout    = 'grid_5x4_lay_square.mat';
    cfg.component = [i];
    cfg.zlim      = 'maxabs';
    cfg.highlight = 'labels';
    
    ft_topoplotIC(cfg, comp)
    colormap('jet')
    
end

%% plot one component 
figure

cfg           = [];
cfg.layout    = 'grid_5x4_lay_square.mat';
cfg.component = 1;
cfg.zlim      = 'maxabs';
cfg.highlight = 'labels';

ft_topoplotIC(cfg, comp)
colormap('jet')

%% compute  
 
%% removing components H2

R2 = h2_remove(Y(1,:),X(1,:));
figure
subplot(size(X,1),1,1);
plot(R2(:,1:2000));
for i = 2 : size(X,1)
    
    R2 = h2_remove(R2,X(i,:));
    subplot(size(X,1),1,i);
    plot(R2(:,1:2000));
end
%% removing linear effect
cfg                = []; 
cfg.bidsFolder     = '/home/matteo/Desktop/tle_e/converted/';
cfg.sitFname       = 'sub-RESP0465_ses-SITUATION1A_task-acute_ieeg';
[nResData,resData] = remove_resected_influence(cfg,m_data);


cfg           = [];
cfg.viewmode  = 'vertical';  
cfg.blocksize = 15;
cfg.seldat    = 'current';

ft_databrowser(cfg,m_data)

ft_databrowser(cfg,nResData)

ft_databrowser(cfg,resData)


%% connectivity analysis

cfg         = [];
cfg.trials  = 'all';
cfg.length  = 5; %seconds of new trials
cfg.overlap = 0;

nResData = ft_redefinetrial(cfg,nResData);
cutData  = ft_redefinetrial(cfg,m_data);

%% 
cfg              = [];
cfg.output       = 'fourier';
cfg.channel      = 'all';
cfg.method       = 'mtmfft';
cfg.taper        = 'hanning';
cfg.foi          = [40:20:150];                         % analysis 2 to 30 Hz in steps of 2 Hz
%cfg.t_ftimwin    = ones(length(cfg.foi),1).*1;     % length of time window = 5 sec
cfg.tapsmofrq    = 2;
cfg.toi          = cellfun(@median,cutData.time);
       
cutTFR           = ft_freqanalysis(cfg, cutData);

nResTFR           = ft_freqanalysis(cfg, nResData);

%%
cfg = [];
cfg.showlabels   = 'yes';
cfg.layout       = 'grid_5x4_lay_square.mat';
figure
ft_topoplotTFR(cfg, cutTFR);
figure
ft_topoplotTFR(cfg, nResTFR);



%%
cfg         = [];
cfg.method  = 'coh';

connVirt     = ft_connectivityanalysis(cfg,nResTFR);

cfg.channel  = nResData.label;
connReal     = ft_connectivityanalysis(cfg,cutTFR);
%%
cfg           = [];
cfg.parameter = 'cohspctrm';
cfg.zlim      = [0 1];
figure
ft_connectivityplot(cfg, connReal,connVirt);
%figure
%ft_connectivityplot(cfg, connVirt);



