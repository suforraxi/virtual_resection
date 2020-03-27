function plot_coh_timecourse()

inFolder = '/home/matteo/Desktop/virtual_resection/syn/h2_all_epochs/';

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
                  if(any(strcmp('syn',fieldnames(sit_res{1}))))
                    coh_tc{s}.pre_syn_v   = cellfun(@getfield,sit_res,repmat({'syn'},1,numel(sit_res)));
                    coh_tc{s}.virt1_syn_v = cellfun(@getfield,sit_res,repmat({'V1syn'},1,numel(sit_res)));
                    coh_tc{s}.virt2_syn_v = cellfun(@getfield,sit_res,repmat({'V2syn'},1,numel(sit_res)));
                  end
                  
                  
                  coh_tc{s}.Cpre        = cellfun(@getfield,sit_res,repmat({'C'},1,numel(sit_res)),'UniformOutput',false);
                   
                  
                  coh_tc{s}.C1           = cellfun(@getfield,sit_res,repmat({'C1'},1,numel(sit_res)),'UniformOutput',false);
                  
                  
                  coh_tc{s}.C2          = cellfun(@getfield,sit_res,repmat({'C2'},1,numel(sit_res)),'UniformOutput',false);
                    
               
            case 'Post'
                 if(any(strcmp('syn',fieldnames(sit_res{1}))))
                    coh_tc{s}.post_syn_v  = cellfun(@getfield,sit_res,repmat({'syn'},1,numel(sit_res)));
                 end
                  coh_tc{s}.Cpost       = cellfun(@getfield,sit_res,repmat({'C'},1,numel(sit_res)),'UniformOutput',false);
        end
        
        
        coh_tc{s}.fname = strcat(sit_res{1}.subjName,...
                                 '_ep_',num2str(sit_res{1}.epoch),...
                                 '_band_',num2str(sit_res{1}.band(1)),'_',num2str(sit_res{1}.band(end)),...
                                 '_fc_',sit_res{1}.fc_type ...
                                 );    
       coh_tc{s}.fc_type = sit_res{1}.fc_type;
    end
end





% plot synchronizability
fig    = zeros(1,numel(coh_tc));

co     = [1 0 0; 0 0 1; 0 0 0; 0 1 0 ; 1 1 0; 0 1 1;1 0 1];
Lwidth = 1.5;


% for i = 1 : numel(coh_tc)
%    
%     fig(i) = figure;
%     set(fig(i),'defaultAxesColorOrder',co)
%     plot(coh_tc{i}.pre_syn_v,'o-','LineWidth',Lwidth)
%     hold
%     %plot(repmat(mean(coh_tc{i}.pre_coh_v),1,length(coh_tc{i}.pre_coh_v)),'r--')
%     plot(coh_tc{i}.post_syn_v,'o-','LineWidth',Lwidth)
%     %plot(repmat(mean(coh_tc{i}.post_coh_v),1,length(coh_tc{i}.post_coh_v)),'b--')
%     plot(coh_tc{i}.virt1_syn_v,'o-','LineWidth',Lwidth)
%     %plot(repmat(mean(coh_tc{i}.virt1_coh_v),1,length(coh_tc{i}.virt1_coh_v)),'k--')
%     plot(coh_tc{i}.virt2_syn_v,'o-','LineWidth',Lwidth)
%     %plot(repmat(mean(coh_tc{i}.virt2_coh_v),1,length(coh_tc{i}.virt2_coh_v)),'g--')
%     
%     
%     ylim([0 1])
%     ylabel('Synchronizability')
%     xlabel('trial')
%     %legend({'pre','pre-avg','post','post-avg','virtNoOrtho','virtNoOrtho-avg','virtual','virtual-avg'})
%     legend({'pre','post','VR_{Naive}','VR_{Ortho}'})
%     
%     set(fig(i),'defaultAxesColorOrder',co)
%     
%     title(subjNames{i})
% end
% % save syncronizability
% outFolder = '/home/matteo/Desktop/virtual_resection/res_pic/syn/';
% if(~exist(outFolder,'dir'))
%     mkdir(outFolder)
% end
% for i = 1 : numel(fig)
% 
%     set(fig(i),'WindowState','fullscreen')
%     outfile = fullfile(outFolder,coh_tc{i}.fname);
%     print(fig(i),outfile,'-djpeg');
%     close(fig(i))
% end

% plot functinal connectivity matrix
% fig    = zeros(1,numel(coh_tc));
% nr     = 4;
% nTrial = numel(coh_tc{1}.Cpre);
%     
% fSizeMatrix = 12;
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
%         max_t(j) = max([max((m_pre(:))),max((m_v1(:))),max((m_v2(:))),max((m_post(:)))]);
%         min_t(j) = min([min((m_pre(:))),min((m_v1(:))),min((m_v2(:))),min((m_post(:)))]);
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
%         imagesc(abs(m_pre),[c_min,c_max]);
%         title(sprintf('Trial %i',j))
%         subplot(nr,nTrial,j+nTrial)
%         imagesc(abs(m_v1),[c_min,c_max]);
%         subplot(nr,nTrial,j+2*nTrial)
%         imagesc(abs(m_v2),[c_min,c_max]);
%         subplot(nr,nTrial,j+3*nTrial)
%         imagesc(abs(m_post),[c_min,c_max]);
%         
%         
%     end
%     
%     subplot(nr,nTrial,1)
%     ylabel(sprintf('%s Pre',subjNames{i}),'FontSize',fSizeMatrix)
%     subplot(nr,nTrial,nTrial+1)
%     ylabel('VR_{naive}','FontSize',fSizeMatrix)
%     subplot(nr,nTrial,2*nTrial+1)
%     ylabel('Vr_{ortho}','FontSize',fSizeMatrix)
%     subplot(nr,nTrial,3*nTrial+1)
%     ylabel('Post','FontSize',fSizeMatrix)
%     
%    
%     colorbar('manual','Position',[0.94 0.1 0.02 0.2])
%    
%     
% end
% 
% 
% % save  functional connectivity matrices
% outFolder = '/home/matteo/Desktop/virtual_resection/res_pic/matrix/';
% if(~exist(outFolder,'dir'))
%     mkdir(outFolder)
% end
% for i = 1 : numel(fig)
% 
%     set(fig(i),'WindowState','fullscreen')
%     outfile = fullfile(outFolder,strcat(coh_tc{i}.fname,'_','matrix'));
%     print(fig(i),outfile,'-djpeg');
%     close(fig(i))
% end


% plot synchronizability
 if(any(strcmp('pre_syn_v',fieldnames(coh_tc{1}))))
    plot_syn_epochs(coh_tc,subjNames,Lwidth,co)
 end
% global connectivity

%plot_global(coh_tc,subjNames,Lwidth,co)
plot_global_all_epochs(coh_tc,subjNames,Lwidth,co)


function plot_syn_epochs(coh_tc,subjNames,Lwidth,co)

%fig     = zeros(1,numel(coh_tc));
nTrial  = numel(coh_tc{1}.Cpre);


nr      = 4;
fig     = figure;
for i = 1 : numel(coh_tc)
   
    
   
    subplot(2,4,i)
    
    pre_v  = coh_tc{i}.pre_syn_v;
    v1_v   = coh_tc{i}.virt1_syn_v;
    v2_v   = coh_tc{i}.virt2_syn_v;
    post_v = coh_tc{i}.post_syn_v;
    
     val    = [pre_v , post_v, v1_v, v2_v];
     labels = [repmat({'pre'},1,numel(pre_v)) , repmat({'post'},1,numel(post_v)), ...
               repmat({'VR_{Naive}'},1,numel(v1_v)) , repmat({'VR_{Ortho}'},1,numel(v2_v)), ...
              ];
    
    violinplot(val,labels);     

    ylabel('Synchronizability')
    title(subjNames{i})
    %ylim([0 1])
    
end

% save global functional connectivity
outFolder = '/home/matteo/Desktop/virtual_resection/res_pic/all_subjects/';
if(~exist(outFolder,'dir'))
    mkdir(outFolder)
end
for i = 1 : numel(fig)

    set(fig(i),'WindowState','fullscreen')
    outfile = fullfile(outFolder,strcat('all_subjects','_',coh_tc{1}.fc_type,'_','syn'));
    print(fig(i),outfile,'-djpeg');
    close(fig(i))
end

function plot_global_all_epochs(coh_tc,subjNames,Lwidth,co)

fig     = zeros(1,numel(coh_tc));
nTrial  = numel(coh_tc{1}.Cpre);


nr      = 4;
fig     = figure;
vari     = zeros(nr,numel(coh_tc));
mu      = zeros(nr,numel(coh_tc));
err     = zeros(nr,numel(coh_tc));
h       = zeros(3,numel(coh_tc));
p       = zeros(3,numel(coh_tc));

relErr = @(m,x,y) m(y,x);%@(m,x,y) abs((m(x)-m(y))/m(y));%abs((mean(x)-mean(y)));
for i = 1 : numel(coh_tc)
   
    
    
    subplot(4,2,i)
    %avg_fc  = zeros(nr,nTrial);
  
    pre_v  = cellfun(@get_Stats,coh_tc{i}.Cpre);
    v1_v   = cellfun(@get_Stats,coh_tc{i}.C1);
    v2_v   = cellfun(@get_Stats,coh_tc{i}.C2);
    post_v = cellfun(@get_Stats,coh_tc{i}.Cpost);
    
    mu(1,i)  = mean(pre_v);
    mu(2,i)  = mean(v1_v);
    mu(3,i)  = mean(v2_v);
    mu(4,i)  = mean(post_v);
    
    
    vari(1,i)  = std(pre_v);
    vari(2,i)  = std(v1_v);
    vari(3,i)  = std(v2_v);
    vari(4,i)  = std(post_v);
    
    
    
    
    [h(1,i),p(1,i)] = kstest2(pre_v,post_v,'Tail','smaller','Alpha',0.01);
    [h(2,i),p(2,i)] = kstest2(v1_v,post_v,'Tail','unequal','Alpha',0.01);
    [h(3,i),p(3,i)] = kstest2(v2_v,post_v,'Tail','unequal','Alpha',0.01);
    
    
    
    val    = [pre_v , post_v, v1_v, v2_v];
    labels = [repmat({'pre'},1,numel(pre_v)) , repmat({'post'},1,numel(post_v)), ...
               repmat({'VR_{Naive}'},1,numel(v1_v)) , repmat({'VR_{partialization}'},1,numel(v2_v)), ...
             ];
    
    violinplot(val,labels,'showMean',true);     
    hold
    line([0 4.5],[0.1 0.1],'LineStyle','--','Color','green') 
    ylabel('Global Connectivity')
    title(subjNames{i})
    %ylim([0 0.45])
    
    
    
    
end

(vari-repmat(vari(1,:),4,1))./repmat(vari(1,:),4,1)*100
(mu-repmat(mu(1,:),4,1))./repmat(mu(1,:),4,1)*100
cv = vari./mu
(cv-repmat(cv(1,:),4,1))./repmat(cv(1,:),4,1)*100


% save global functional connectivity
outFolder = '/home/matteo/Desktop/virtual_resection/res_pic/all_epochs/';
if(~exist(outFolder,'dir'))
    mkdir(outFolder)
end
if(numel(fig) == 1)
        set(fig,'WindowState','fullscreen')
        outfile = fullfile(outFolder,strcat('all_subjects','_',coh_tc{1}.fc_type,'_','global'));
        print(fig,outfile,'-djpeg');
        close(fig)
else
    for i = 1 : numel(fig)

        set(fig(i),'WindowState','fullscreen')
        outfile = fullfile(outFolder,strcat(coh_tc{i}.fname,'_','global'));
        print(fig(i),outfile,'-djpeg');
        close(fig(i))
    end
end

function plot_global(coh_tc,subjNames,Lwidth,co)

fig     = zeros(1,numel(coh_tc));
nTrial  = numel(coh_tc{1}.Cpre);


nr      = 4;
for i = 1 : numel(coh_tc)
   
    fig(i) = figure;
    subplot(1,2,1)
    %max_t = zeros(1,nTrial);
    %min_t = zeros(1,nTrial);
    
    avg_fc  = zeros(nr,nTrial);
  
    
    for j = 1 : nTrial
        
        m_pre  = coh_tc{i}.Cpre{j};
        m_v1   = coh_tc{i}.C1{j};
        m_v2   = coh_tc{i}.C2{j};
        m_post = coh_tc{i}.Cpost{j};
       
        avg_fc(1,j) = get_Stats(m_pre);
        avg_fc(2,j) = get_Stats(m_post);
        avg_fc(3,j) = get_Stats(m_v1);
        avg_fc(4,j) = get_Stats(m_v2);
        
    end
    
    set(fig(i),'defaultAxesColorOrder',co)
    plot(avg_fc','o-','LineWidth',Lwidth);
    legend({'pre','post','VR_{Naive}','VR_{Ortho}'});
    title(subjNames{i})
    ylabel('Global Connectivity')
    xlabel('Trials')
    ylim([0 1])
  
    subplot(1,2,2)
    notBoxPlot(avg_fc')
    xticklabels({'pre','post','VR_{Naive}','VR_{Ortho}'});
    ylabel('Global Connectivity')
    ylim([0 0.4])
    
end

% save global functional connectivity
outFolder = '/home/matteo/Desktop/virtual_resection/res_pic/global/';
if(~exist(outFolder,'dir'))
    mkdir(outFolder)
end
for i = 1 : numel(fig)

    set(fig(i),'WindowState','fullscreen')
    outfile = fullfile(outFolder,strcat(coh_tc{i}.fname,'_','global'));
    print(fig(i),outfile,'-djpeg');
    close(fig(i))
end


function [avg_fc]= get_Stats(m)

        % compute average global connectivity
        get_Gfc = @(x) sum(sum(x))/(size(x,1)^2-size(x,1));
        
        % put zeros on the diagonal
        m(1:length(m)+1:end)   = 0;
        
        avg_fc = get_Gfc(abs(m)); % FIX abs value if it is correlation 
        
        %mask = true(size(m));
        
        %mask(1:length(m)+1:end) = false;
        
        %var_fc = std(m(mask));
       
        
        
        
        