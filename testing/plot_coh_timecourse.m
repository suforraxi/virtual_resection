
inFolder = '/home/matteo/Desktop/virtual_resection/coh/';

inFiles  = dir(fullfile(inFolder,'*.mat'));

subjNames = unique(extractBetween({inFiles(:).name},5,12));


getCoh   = @(x) x.coh;
getVCoh  = @(x) x.Vcoh;

coh_tc   = cell(1,numel(subjNames)); 
for s = 1: numel(subjNames)
    
    subjFiles      = dir(fullfile(inFolder,strcat('*',subjNames{s},'*')));
  
    for f = 1 : numel(subjFiles)
        
        load(fullfile(subjFiles(f).folder,subjFiles(f).name));
        
        %coh_tc{s}.fname = sit_name
        switch sit_res{1}.sitType
           
            case 'Pre'
                  coh_tc{s}.pre_coh_v  = cellfun(getCoh,sit_res);
                   
                  coh_tc{s}.virt_coh_v = cellfun(getVCoh,sit_res);
                    
               
            case 'Post'
                   coh_tc{s}.post_coh_v = cellfun(getCoh,sit_res);
            
        end
        coh_tc{s}.fname = strcat(sit_res{1}.subjName,...
                                 '_ep_',num2str(sit_res{1}.epoch),...
                                 '_band_',num2str(sit_res{1}.band(1)),'_',num2str(sit_res{1}.band(end)) ...
                                 );     
    end
end

fig = zeros(1,numel(coh_tc));
for i = 1 : numel(coh_tc)
   
    fig(i) = figure;
    plot(coh_tc{i}.pre_coh_v,'ro-')
    hold
    plot(repmat(mean(coh_tc{i}.pre_coh_v),1,length(coh_tc{i}.pre_coh_v)),'r--')
    plot(coh_tc{i}.post_coh_v,'bo-')
    plot(repmat(mean(coh_tc{i}.post_coh_v),1,length(coh_tc{i}.post_coh_v)),'b--')
    plot(coh_tc{i}.virt_coh_v,'go-')
    plot(repmat(mean(coh_tc{i}.virt_coh_v),1,length(coh_tc{i}.virt_coh_v)),'g--')
    ylim([0 1])
    ylabel('Coherence')
    xlabel('trial')
    legend({'pre','pre-avg','post','post-avg','virtual','virtual-avg'})
    title(subjNames{i})
end

outFolder = '/home/matteo/Desktop/virtual_resection/res_pic/';
for i = 1 : numel(fig)

    set(fig(i),'WindowState','fullscreen')
    outfile = fullfile(outFolder,coh_tc{i}.fname);
    print(fig(i),outfile,'-djpeg');
    close(fig(i))
end


