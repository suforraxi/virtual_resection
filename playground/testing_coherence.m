path settings

addpath('/home/matteo/Desktop/git_rep/virtual_resection/')
addpath('/home/matteo/Desktop/git_rep/virtual_resection/layout_utils/')
addpath('/home/matteo/Desktop/git_rep/virtual_resection/testing/')


modulation - message & carrier

fs = 1000;              % sampling frequency
Ts = 1/fs;              % sampling period 
N  = 1000;              % number of samples
t  = [0:Ts:N*Ts-Ts];    % time axis

fc = 120;               % carrier frequency
fm = 12;                % modulation frequency
p  = sin(2*pi*fc*t);    % carrier signal
x  = sin(2*pi*fm*t);    % modulation signal
m1  = 0.75;              % modulation index
m2  = 0.45;

y1  = ((1 + m1*x).*p);      % modulated signal
y2  = ((1 + m2*x).*p);      % modulated signal

figure
subplot(411)
plot(p)
title('Carrier')
subplot(412)
plot(x)
title('Message')
subplot(413)
plot(y1)
title(sprintf('modulation %1.2f',m1))
subplot(414)
plot(y2)
title(sprintf('modulation %1.2f',m2))




create example data
nch  = 31;
N    = length(y);

data = repmat(y,nch,1);

a          = [-0.6] ;
b          = [0.6 ];

c          =  0:0.2:6;

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
view data
cfg = []; cfg.viewmode = 'vertical'; cfg.blocksize = 15; cfg.seldat = 'current';
ft_databrowser(cfg,m_data)
plot components
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

