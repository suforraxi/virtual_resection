
inFolder = '/home/matteo/Desktop/virtual_resection/coh/';

inFiles  = dir(fullfile(inFolder,'*.mat'));

subjNames = unique(extractBetween({inFiles(:).name},5,12));

%getCoh   = @(x) x.coh;
%getV1Coh  = @(x) x.V1coh;
%getV2Coh  = @(x) x.V2coh;

coh_tc   = cell(1,numel(subjNames)); 
for s = 1: numel(subjNames)
    
    subjFiles      = dir(fullfile(inFolder,strcat('*',subjNames{s},'*')));
  
    for f = 1 : numel(subjFiles)
        
        load(fullfile(subjFiles(f).folder,subjFiles(f).name));
        
        %coh_tc{s}.fname = sit_name
        switch sit_res{1}.sitType
           
            case 'Pre'
                  coh_tc{s}.pre_coh_v   = cellfun(@getfield,sit_res,repmat({'coh'},1,numel(sit_res)));
                  coh_tc{s}.Cpre        = cellfun(@getfield,sit_res,repmat({'C'},1,numel(sit_res)),'UniformOutput',false);
                   
                  coh_tc{s}.virt1_coh_v = cellfun(@getfield,sit_res,repmat({'V1coh'},1,numel(sit_res)));
                  coh_tc{s}.C1           = cellfun(@getfield,sit_res,repmat({'C1'},1,numel(sit_res)),'UniformOutput',false);
                  
                  coh_tc{s}.virt2_coh_v = cellfun(@getfield,sit_res,repmat({'V2coh'},1,numel(sit_res)));
                  coh_tc{s}.C2          = cellfun(@getfield,sit_res,repmat({'C2'},1,numel(sit_res)),'UniformOutput',false);
                    
               
            case 'Post'
                  coh_tc{s}.post_coh_v  = cellfun(@getfield,sit_res,repmat({'coh'},1,numel(sit_res)));
                  coh_tc{s}.Cpost       = cellfun(@getfield,sit_res,repmat({'C'},1,numel(sit_res)),'UniformOutput',false);
        end
        coh_tc{s}.fname = strcat(sit_res{1}.subjName,...
                                 '_ep_',num2str(sit_res{1}.epoch),...
                                 '_band_',num2str(sit_res{1}.band(1)),'_',num2str(sit_res{1}.band(end)),...
                                 '_fc_',sit_res{1}.fc_type ...
                                 );     
    end
end

% plot synchronizability
fig    = zeros(1,numel(coh_tc));

co     = [1 0 0; 0 0 1; 0 0 0; 0 1 0 ; 1 1 0; 0 1 1;1 0 1];
Lwidth = 1.5;
for i = 1 : numel(coh_tc)
   
    fig(i) = figure;
    set(fig(i),'defaultAxesColorOrder',co)
    plot(coh_tc{i}.pre_coh_v,'o-','LineWidth',Lwidth)
    hold
    %plot(repmat(mean(coh_tc{i}.pre_coh_v),1,length(coh_tc{i}.pre_coh_v)),'r--')
    plot(coh_tc{i}.post_coh_v,'o-','LineWidth',Lwidth)
    %plot(repmat(mean(coh_tc{i}.post_coh_v),1,length(coh_tc{i}.post_coh_v)),'b--')
    plot(coh_tc{i}.virt1_coh_v,'o-','LineWidth',Lwidth)
    %plot(repmat(mean(coh_tc{i}.virt1_coh_v),1,length(coh_tc{i}.virt1_coh_v)),'k--')
    plot(coh_tc{i}.virt2_coh_v,'o-','LineWidth',Lwidth)
    %plot(repmat(mean(coh_tc{i}.virt2_coh_v),1,length(coh_tc{i}.virt2_coh_v)),'g--')
    
    
    ylim([0 1])
    ylabel('Synchronizability')
    xlabel('trial')
    %legend({'pre','pre-avg','post','post-avg','virtNoOrtho','virtNoOrtho-avg','virtual','virtual-avg'})
    legend({'pre','post','VR_{Naive}','VR_{Ortho}'})
    
    set(fig(i),'defaultAxesColorOrder',co)
    
    title(subjNames{i})
end
% save syncronizability
outFolder = '/home/matteo/Desktop/virtual_resection/res_pic/';
for i = 1 : numel(fig)

    set(fig(i),'WindowState','fullscreen')
    outfile = fullfile(outFolder,coh_tc{i}.fname);
    print(fig(i),outfile,'-djpeg');
    close(fig(i))
end

% plot functinal connectivity matrix
% fig    = zeros(1,numel(coh_tc));
% nr     = 4;
% nTrial = numel(coh_tc{1}.Cpre);
%     
% 
% 
% for i = 1 : numel(coh_tc)
%    
%     fig(i) = figure;
%     max_t = zeros(1,nTrial);
%     min_t = zeros(1,nTrial);
%     
%     for j = 1 : nTrial
%         
%         m_pre  = coh_tc{i}.Cpre{j};
%         m_v1   = coh_tc{i}.C1{j};
%         m_v2   = coh_tc{i}.C2{j};
%         m_post = coh_tc{i}.Cpost{j};
%         
%         m_pre(1:length(m_pre)+1:end)   = 0;
%         m_v1(1:length(m_v1)+1:end)     = 0;
%         m_v2(1:length(m_v2)+1:end)     = 0;
%         m_post(1:length(m_post)+1:end) = 0;
%         
%         max_t(j) = max([max(m_pre(:)),max(m_v1(:)),max(m_v2(:)),max(m_post(:))]);
%         min_t(j) = min([min(m_pre(:)),min(m_v1(:)),min(m_v2(:)),min(m_post(:))]);
%     end
%     
%     c_max = max(max_t); 
%     c_min = min(min_t);
%     
%     for j = 1 : nTrial
%         
%         m_pre  = coh_tc{i}.Cpre{j};
%         m_v1   = coh_tc{i}.C1{j};
%         m_v2   = coh_tc{i}.C2{j};
%         m_post = coh_tc{i}.Cpost{j};
%         
%         m_pre(1:length(m_pre)+1:end)   = 0;
%         m_v1(1:length(m_v1)+1:end)     = 0;
%         m_v2(1:length(m_v2)+1:end)     = 0;
%         m_post(1:length(m_post)+1:end) = 0;
%         
%         subplot(nr,nTrial,j)
%         imagesc(m_pre,[c_min,c_max]);
%         subplot(nr,nTrial,j+nTrial)
%         imagesc(m_v1,[c_min,c_max]);
%         subplot(nr,nTrial,j+2*nTrial)
%         imagesc(m_v2),[c_min,c_max];
%         subplot(nr,nTrial,j+3*nTrial)
%         imagesc(m_post,[c_min,c_max]);
%         
%         
%     end
%     
% end

% global connectivity

% plot functional connectivity matrix
fig     = zeros(1,numel(coh_tc));
nTrial  = numel(coh_tc{1}.Cpre);
get_Gfc = @(x) sum(sum(x))/(size(x,1)^2-size(x,1));    
nr      = 4;
for i = 1 : numel(coh_tc)
   
    fig(i) = figure;
    max_t = zeros(1,nTrial);
    min_t = zeros(1,nTrial);
    
    g_fc  = zeros(nr,nTrial);
    for j = 1 : nTrial
        
        m_pre  = coh_tc{i}.Cpre{j};
        m_v1   = coh_tc{i}.C1{j};
        m_v2   = coh_tc{i}.C2{j};
        m_post = coh_tc{i}.Cpost{j};
        
        m_pre(1:length(m_pre)+1:end)   = 0;
        m_v1(1:length(m_v1)+1:end)     = 0;
        m_v2(1:length(m_v2)+1:end)     = 0;
        m_post(1:length(m_post)+1:end) = 0;
        
        g_fc(1,j) = get_Gfc(abs(m_pre));
        g_fc(2,j) = get_Gfc(abs(m_post));
        g_fc(3,j) = get_Gfc(abs(m_v1));
        g_fc(4,j) = get_Gfc(abs(m_v2));
        
    end
    
    set(fig(i),'defaultAxesColorOrder',co)
    plot(g_fc','o-','LineWidth',Lwidth);
    legend({'pre','post','VR_{Naive}','VR_{Ortho}',});
    title(subjNames{i})
    ylabel('Global Connectivity')
    xlabel('Trials')
    ylim([0 1])
   
end

% save global functional connectivity

for i = 1 : numel(fig)

    set(fig(i),'WindowState','fullscreen')
    outfile = fullfile(outFolder,strcat(coh_tc{i}.fname,'_','global'));
    print(fig(i),outfile,'-djpeg');
    close(fig(i))
end

