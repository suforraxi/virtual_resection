function testing_coherence()
addpath('/home/matteo/Desktop/git_rep/virtual_resection/')
close all

fs = 1000;              % sampling frequency
Ts = 1/fs;              % sampling period 
N  = 1000;              % number of samples
t  = [0:Ts:N*Ts-Ts];    % time axis

fc = 120;               % carrier frequency
fm = 12;                % modulation frequency
p  = sin(2*pi*fc*t);    % carrier signal
x  = sin(2*pi*fm*t);    % modulation signal
m  = 0.75;              % modulation index

y  = ((1 + m*x).*p);      % modulated signal

nch  = 31;
N    = length(y);

data = repmat(y,nch,1);

a          = [-0.6] ;
b          = [0.6 ];

c          =  0:0.2:6;



% create data
d = [];
load('../layout_utils/m_data.mat');

for j = 1 : numel(a)
    
    d{j}.label   = m_data.label(1:nch); 
    d{j}.fsample = fs;
    d{j}.time    = {1:N};

    for i = 1 : nch 
        noise = (a(j) + (b(j)-a(j)).*rand(1,N));
        data(i,:) = data(i,:) + noise*c(i);
    end

    %d{j}.trial = {get_Envelope(data)};
    d{j}.trial = {data};
end





% plot data

for i = 1 : numel(d)
   
    data = d{i}.trial{1};
    nch  = size(data,1);
    dEnv = get_Envelope(data);
    figure
    for j = 1 : nch
       
        subplot(nch,1,j);
        plot(data(j,:),'b')
        %hold
        %plot(dEnv(j,:),'g')
    end

end

%  cfg           = [];
%  cfg.viewmode  = 'vertical';  
%  cfg.blocksize = 15;
%  cfg.seldat    = 'current';
%  
%  ft_databrowser(cfg,d{1})
% 



%% ICA or PCA
% cfg              = [];
% cfg.method       = 'pca';
% cfg.channel      = 'all';
% cfg.trials       = 'all';
% cfg.numcomponent = 'all';
% cfg.demean       = 'yes';
% cfg.updatesens   = 'yes';
% cfg.feedback     = 'text';
% 
% 
% comp = ft_componentanalysis(cfg,d{1});

%% browse ICA

 %cfg           = [];
 %cfg.layout    = 'grid_5x4_lay_square.mat';
 %cfg.zlim      = 'maxabs';
 %cfg.highlight = 'labels';
 %ft_icabrowser(cfg,comp)

%   figure
%  for i = 1 : numel(comp.label)
%      
%      subplot(8,4,i)
%      
%      cfg           = [];
%      cfg.layout    = 'grid_5x4_lay_square.mat';
%      cfg.component = [i];
%      cfg.zlim      = 'maxabs';
%      cfg.style     = 'straight';
%      %cfg.highlight = 'labels';
%      cfg.interplimits = 'electrodes';
%      
%      ft_topoplotIC(cfg, comp)
%      colormap('jet')
%      
%  end
%  
 
 % compute correlation and coherence

for i =  1 : 1%numel(d)
   
    cfg              = [];
    cfg.demean       = 'yes';
    d{i}             = ft_preprocessing(cfg,d{i});
    data = d{i}.trial{1};
    nch  = size(data,1);
   
    
    mEnv            = get_Envelope(m_data.trial{i});
    C               = fc(mEnv);
    [coh,C,V,E]     = get_coherence(C);
    
    plot_AC(d{i},coh,C,V);


end

vr              = remove_ch(d{1},1);
[vcoh,vC,vV,vE] = get_coherence(vr.trial{1});
plot_AC(vr,vcoh,vC,vV)

vr;
% function plot_AC(d,coh,C,V)
%     
%     fig_C = figure;
%     fig_V = figure;
%     figure(fig_C)
%     C(1:size(C,1)+1:end) = 0; 
%     imagesc(C)
%     title(sprintf('coh:%.4f',coh))
%     
% 
%     figure(fig_V)
%     plot_topo(d,V);
% 
% 
% function d = remove_ch(d,idx_ch)
% 
%     nTrial      = numel(d.trial); 
%     idx         = true(length(d.label),1);
%     idx(idx_ch) = false;
%     
%     for i = 1 : nTrial 
%         m          = d.trial{i};
%         m          = get_ortho_matrix(m,idx_ch);
%         d.trial{i} = m;
%         d.label    = d.label(idx);   
%         
%     end
% 
% 
% 
%  % topoplot components
% function plot_topo(d,V) 
%   
%     mycomp           = d;
%     mycomp.topo      = V;
%     mycomp.unmixing  = inv(V);
%     mycomp.topolabel = d.label;
%     
%     nr = ceil(sqrt(length(V)));
%     nc = nr;
%     for j = 1 : numel(mycomp.label)
%      
%      subplot(nr,nc,j)
%      
%      cfg           = [];
%      cfg.layout    = 'grid_5x4_lay_square.mat';
%      cfg.component = [j];
%      cfg.zlim      = 'maxabs';
%      %cfg.highlight = 'labels';
%      cfg.interplimits = 'electrodes';
%      cfg.style     = 'straight';
%    
%      ft_topoplotIC(cfg, mycomp)
%      colormap('jet')
%      
%     end
     

%% view data
% cfg           = [];
% cfg.viewmode  = 'vertical';  
% cfg.blocksize = 15;
% cfg.seldat    = 'current';
% 
% ft_databrowser(cfg,m_data)

%cfg=[];
%cfg.refchannel={'all'}; 
%m_data = ft_preprocessing(cfg,m_data);

% %% ICA or PCA
% cfg              = [];
% cfg.method       = 'fastica';
% cfg.channel      = 'all';
% cfg.trials       = 'all';
% cfg.numcomponent = 'all';
% cfg.demean       = 'yes';
% cfg.updatesens   = 'yes';
% cfg.feedback     = 'text';
% 
% 
% comp = ft_componentanalysis(cfg,m_data);
% 
%% plot components



%  

% %% plot compenent and time-course
% cfg            = [];
% cfg.layout     = 'grid_5x4_lay_square.mat';
% cfg.viewmode   = 'component';
% cfg.continuous = 'yes';
% cfg.blocksize  = 60;
% cfg.channels   = 'all';
% cfg.zlim       = 'maxabs';
% ft_databrowser(cfg,comp);

%% topoplot of components
% figure
% for i = 1 : numel(comp.label)
%     
%     subplot(5,4,i)
%     
%     cfg           = [];
%     cfg.layout    = 'grid_5x4_lay_square.mat';
%     cfg.component = [i];
%     cfg.zlim      = 'maxabs';
%     cfg.highlight = 'labels';
%     
%     ft_topoplotIC(cfg, comp)
%     colormap('jet')
%     
% end

% %% plot one component 
% figure
% 
% cfg           = [];
% cfg.layout    = 'grid_5x4_lay_square.mat';
% cfg.component = 1;
% cfg.zlim      = 'maxabs';
% cfg.highlight = 'labels';
% 
% ft_topoplotIC(cfg, comp)
% colormap('jet')

%% compute  
 


% cfg           = [];
% cfg.viewmode  = 'vertical';  
% cfg.blocksize = 15;
% cfg.seldat    = 'current';
% 
% ft_databrowser(cfg,m_data)
% 
% ft_databrowser(cfg,nResData)
% 
% ft_databrowser(cfg,resData)


%% connectivity analysis

% cfg         = [];
% cfg.trials  = 'all';
% cfg.length  = 5; %seconds of new trials
% cfg.overlap = 0;
% 
% nResData = ft_redefinetrial(cfg,nResData);
% cutData  = ft_redefinetrial(cfg,m_data);
% 
% %% 
% cfg              = [];
% cfg.output       = 'fourier';
% cfg.channel      = 'all';
% cfg.method       = 'mtmfft';
% cfg.taper        = 'hanning';
% cfg.foi          = [40:20:150];                         % analysis 2 to 30 Hz in steps of 2 Hz
% %cfg.t_ftimwin    = ones(length(cfg.foi),1).*1;     % length of time window = 5 sec
% cfg.tapsmofrq    = 2;
% cfg.toi          = cellfun(@median,cutData.time);
%        
% cutTFR           = ft_freqanalysis(cfg, cutData);
% 
% nResTFR           = ft_freqanalysis(cfg, nResData);

%%
% cfg = [];
% cfg.showlabels   = 'yes';
% cfg.layout       = 'grid_5x4_lay_square.mat';
% figure
% ft_topoplotTFR(cfg, cutTFR);
% figure
% ft_topoplotTFR(cfg, nResTFR);
% 


%%
% cfg         = [];
% cfg.method  = 'coh';
% 
% connVirt     = ft_connectivityanalysis(cfg,nResTFR);
% 
% cfg.channel  = nResData.label;
% connReal     = ft_connectivityanalysis(cfg,cutTFR);
% %%
% cfg           = [];
% cfg.parameter = 'cohspctrm';
% cfg.zlim      = [0 1];
% figure
% ft_connectivityplot(cfg, connReal,connVirt);
%figure
%ft_connectivityplot(cfg, connVirt);





