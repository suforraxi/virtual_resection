% Save figures and tables of the analysis
%
%  see Demuru et al. doi: 
%

%     Copyright (C) 2020 Matteo Demuru
% 
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <https://www.gnu.org/licenses/>.

function results_and_figures()

% input folder of the results (see virtual_resection_main)
inFolder  = './output/syn/h2/';
% output folder where to save the figures and tables
outFolder = './output/';

% figure Name 
figFname = 'all_epochs';
resFname = 'summary_results.txt';
tabFname = 'kstest_pval_tlb.txt';


inFiles  = dir(fullfile(inFolder,'*.mat'));

subjNames = unique(extractBetween({inFiles(:).name},5,12));


fc_res   = cell(1,numel(subjNames)); 
for s = 1: numel(subjNames)
    
    subjFiles      = dir(fullfile(inFolder,strcat('*',subjNames{s},'*')));
  
    for f = 1 : numel(subjFiles)
        
        load(fullfile(subjFiles(f).folder,subjFiles(f).name));
        
       
        switch sit_res{1}.sitType
           
            case 'Pre'
                  
                  fc_res{s}.Cpre        = cellfun(@getfield,sit_res,repmat({'C'},1,numel(sit_res)),'UniformOutput',false);                 
                  fc_res{s}.C1          = cellfun(@getfield,sit_res,repmat({'C1'},1,numel(sit_res)),'UniformOutput',false);
                  fc_res{s}.C2          = cellfun(@getfield,sit_res,repmat({'C2'},1,numel(sit_res)),'UniformOutput',false);
                    
            case 'Post'
                 
                  fc_res{s}.Cpost       = cellfun(@getfield,sit_res,repmat({'C'},1,numel(sit_res)),'UniformOutput',false);
                  fc_res{s}.so          = sit_res{1}.so;
        end
        
        
        fc_res{s}.fname = strcat(sit_res{1}.subjName,...
                                 '_ep_',num2str(sit_res{1}.epoch),...
                                 '_band_',num2str(sit_res{1}.band(1)),'_',num2str(sit_res{1}.band(end)),...
                                 '_fc_',sit_res{1}.fc_type ...
                                 );    
       fc_res{s}.fc_type = sit_res{1}.fc_type;
    end
end



cfg.fc_res        = fc_res;
cfg.subjNames     = subjNames;
cfg.outFolder     = outFolder;
cfg.thresoldType  = 'prc';
cfg.prctile       = 95;
cfg.alpha_level   = 0.05;
cfg.tabFname      = 'result_pval.tsv';
cfg.resFname      = 'summary_res.tsv';
%plot_global_all_epochs(fc_res,subjNames,outFolder,figFname,resFname,tabFname)
plot_global_fc(cfg)

%  plot and save results:
% 
% 1) Main figure (figFname) showing fot each subject the comparision between the distribuion of global functional connectivity
%    for the four condition: 
%      - pre-resection
%      - post-resection
%      - virtual resection after partialization
%      - naive virtual resection 
%   
% 2) Text file (resFname) summarizing: 
%     - the number of subjects for whom it was not possible formulate a prediction
%     - number of successful predictions for naive virtual resection
%     - number of successful predictions for virtual resection after
%       partialization
%       
% 3) text file (tabFname) witha table containing:
%                - as variable names patients names for whom it was possible to make a prediction
%                - as rows the three tested conditions for the patients for whom it was possible to formulate a prediction
%                  (pre-resection vs post-resection, pre-resection vs naive virtual, pre-resection vs virtual with partialization )
%                - as values the p-values of Kolmogorov-Smirnov test
%
%
%
%
% INPUT
%       fc_res    cell with the functional connectivity result struct for each
%                 patient. Each functional connectivity result struct
%                 contains the following fields:
%             
%                 fname:   string with the coded patient name 
%                          example 'RESPXXXX_ep_(epochlength)_band_(highpass)_(lowpass)_fc_(functionalConnectivityName)'
%                 fc_type: string with the functional connectivity name
%                 Cpre:    cell with the value of global functional
%                          connectivity per epoch during pre-resection recordings
%                 C1:      cell with the value of global functional
%                          connectivity per epoch after naive virtual resection
%                 C2:      cell with the value of global functional
%                          connectivity per epoch after virtual resection with
%                          partialization
%                 Cpost:   cell with the value of global functional
%                          connectivity per epoch during post-resection recordings
%                 so:      seizure outcome string {i.e. EngelClass_AED_InfoAboutMedication:'1A_AED_stop' }               
%       
%       subjNames cell array with the coded names of the patients   
%       outFolder string with the output folder where to save the figures
%                 and tables
%       figFname string, name of the file for the main figure
%       resFname string, name of the file with results in terms of allowed
%                prediction, successful predition for the two type of virtual
%                resections
%       tabFname string, name of the file where the p-value table will be
%                saved
% SAVED FILE
%
%       
%       figFname string, name of the file for the main figure
%       resFname string, name of the file with results in terms of allowed
%                prediction, successful predition for the two type of virtual
%                resections
%       tabFname string, name of the file where the p-value table will be
%                saved

function plot_global_all_epochs(fc_res,subjNames,outFolder,figFname,resFname,tabFname)

nSubjs   = numel(fc_res);

nr      = 4;
nc      = 3;

fig     = figure;
X_lim   = [0 5];
Y_lim   = [0 0.45];

alpha_level = 0.01; % alpha level for Kolmogorov-Sminrnov test

vari    = zeros(nr,numel(fc_res));
mu      = zeros(nr,numel(fc_res));


% matrixes to keep track of the Kolmogorov-Smirnov test
h       = zeros(3,numel(fc_res));
p       = zeros(3,numel(fc_res));



% compute reference threshold for global connectivity magnitude and global
% connectivity variation 
M   = zeros(1,nSubjs); % Magnitude
V   = zeros(1,nSubjs); % Variation
so  = cell(1,nSubjs);
for s = 1 : nSubjs
  
    post_v = cellfun(@get_Stats,fc_res{s}.Cpost);
    
    M(s)  = max(post_v);
    V(s)  = std(post_v);
    so(s) = fc_res{s}.so;
end

% indexes of bad seizure outcome patients (i.e. > Engel 1A )
idx_bad_out = cellfun(@isempty,regexp(so,'1A\w*'));

ref_mag  = max(M(~idx_bad_out));
ref_vari = max(V(~idx_bad_out));

% indexes of subject for whom it is possible to formulate a prediction
answer_allowed = zeros(1,nSubjs);
% indexes of subject for whom naive virtual resection can predict the outcome
virtual1_works = zeros(1,nSubjs);
% indexes of subject for whom virtual resection can predict the outcome
virtual2_works = zeros(1,nSubjs);

for s = 1 : nSubjs
   
   
    subplot(nr,nc,s)
  
    pre_v  = cellfun(@get_Stats,fc_res{s}.Cpre); % pre-resection
    v1_v   = cellfun(@get_Stats,fc_res{s}.C1);   % naive virtual resection
    v2_v   = cellfun(@get_Stats,fc_res{s}.C2);   % virtual resection with partialization
    post_v = cellfun(@get_Stats,fc_res{s}.Cpost); % post-resection
    
    % global connectivity of for different conditions 
    % pre-resection
    mu(1,s)  = mean(pre_v);
    % naive virtual resection
    mu(2,s)  = mean(v1_v);
    % virtual resection with partialization
    mu(3,s)  = mean(v2_v);
    % post-resection
    mu(4,s)  = mean(post_v);
    
    % standard deviation of for different conditions 
    % pre-resection
    vari(1,s)  = std(pre_v);
    % naive virtual resection
    vari(2,s)  = std(v1_v);
    % virtual resection with partialization
    vari(3,s)  = std(v2_v);
    % post-resection
    vari(4,s)  = std(post_v);
    
    
    % individual Kolmogorov-Smirnov test comparsion between different
    % condition:
    %
    % pre-resection distribution vs post-resection distribution
    [h(1,s),p(1,s)] = kstest2(pre_v,post_v,'Tail','smaller','Alpha',alpha_level);
    % pre-resection distribution vs naive virtual resection distribution
    [h(2,s),p(2,s)] = kstest2(pre_v,v1_v,'Tail','smaller','Alpha',alpha_level);
    % pre-resection distribution vs virtual resection after partialization distribution
    [h(3,s),p(3,s)] = kstest2(pre_v,v2_v,'Tail','smaller','Alpha',alpha_level);
    
    
    
    val    = [pre_v , post_v, v1_v, v2_v];
    labels = [repmat({'pre'},1,numel(pre_v)) , repmat({'post'},1,numel(post_v)), ...
               repmat({'VR_{Naive}'},1,numel(v1_v)) , repmat({'VR_{partialization}'},1,numel(v2_v)), ...
             ];
    
    violinplot(val,labels,'showMean',true);     
    
    hold
    
    ref_x  = 0:numel(unique(labels))+1;
    ref_y  = repmat(ref_mag,1,numel(unique(labels))+2);
    errBar = repmat(ref_vari,1,numel(unique(labels))+2);
    
    shadedErrorBar(ref_x, ref_y,errBar, 'lineprops', '-g')
    
    line([0 numel(unique(labels))+1],[ref_mag ref_mag],'LineStyle','--','color','k')
    
    ylabel('Global Connectivity')
    title(subjNames{s})
    
    ylim(Y_lim)
    set(gca,'XLim',X_lim)
    
    if(any( pre_v > ref_mag+ref_vari )) % any value bigger than the reference interval
        answer_allowed(s) = 1;
        
        if(h(2,s))
            virtual1_works(s) = 1;
        end
        if(h(3,s))
            virtual2_works(s) = 1;
        end
    end
end


% save global functional connectivity figures and tables
h = h(:,logical(answer_allowed));
p = p(:,logical(answer_allowed));

% save kstest p-value table
res_tbl = array2table(round(p,5),'VariableNames',subjNames(logical(answer_allowed)),'RowNames',{'pre_vs_post','VRnaive_vs_pre','VRpartialization_vs_pre'});

writetable(res_tbl,fullfile(outFolder,tabFname),'FileType','text','Delimiter','tab','WriteRowNames',1);

fid = fopen(fullfile(outFolder,resFname),'w');
  
fprintf(fid,'successful prediction naive: %i / %i\n',sum(virtual1_works),sum(answer_allowed));
fprintf(fid,'successful prediction partialization: %i / %i\n',sum(virtual2_works),sum(answer_allowed));
fprintf(fid,'cannot answer : %i / %i',sum(~answer_allowed),numel(answer_allowed));

fclose(fid);


if(~exist(outFolder,'dir'))
    mkdir(outFolder)
end
if(numel(fig) == 1)
        set(fig,'WindowState','fullscreen')
        outfile = fullfile(outFolder,strcat(figFname,'_',fc_res{1}.fc_type,'_','global'));
        print(fig,outfile,'-djpeg');
        close(fig)
else
    for s = 1 : numel(fig)

        set(fig(s),'WindowState','fullscreen')
        outfile = fullfile(outFolder,strcat(fc_res{s}.fname,'_','global'));
        print(fig(s),outfile,'-djpeg');
        close(fig(s))
    end
end




function plot_global_fc(cfg)


fc_res      = cfg.fc_res;
subjNames   = cfg.subjNames;
outFolder   = cfg.outFolder;
alpha_level = cfg.alpha_level; % alpha level for Kolmogorov-Sminrnov test
tabFname    = cfg.tabFname;
resFname    = cfg.resFname;
nSubjs   = numel(fc_res);


% compute reference threshold for global connectivity magnitude and global
% connectivity variation 
M   = zeros(1,nSubjs); % Magnitude
V   = zeros(1,nSubjs); % Variation
so  = cellfun(@(x) x.so,fc_res);

% indexes of bad seizure outcome patients (i.e. > Engel 1A )
idx_bad_out = cellfun(@isempty,regexp(so,'1A\w*'));

ref_mag  = 0;
ref_vari = 0;

switch cfg.thresoldType 
    case 'max'
        for s = 1 : nSubjs
  
            post_v = cellfun(@get_Stats,fc_res{s}.Cpost);
    
            M(s)  = max(post_v);
            V(s)  = std(post_v);
            so(s) = fc_res{s}.so;
        end

        ref_mag  = max(M(~idx_bad_out));
        ref_vari = max(V(~idx_bad_out));
    case 'prc'
        post_v = [];
        
        for s = 1 : nSubjs
            if(contains(fc_res{s}.so,'1A_AED'))
                post_v = [ post_v cellfun(@get_Stats,fc_res{s}.Cpost)];
            end
        end
        
        ref_mag = prctile(post_v,cfg.prctile);
end


% indexes of subject for whom it is possible to formulate a prediction
answer_allowed   = zeros(1,nSubjs);
estimation_fc_v1 = zeros(1,nSubjs);
estimation_fc_v2 = zeros(1,nSubjs);

% matrixes to keep track of the Kolmogorov-Smirnov test
h       = -1*ones(8,numel(fc_res));
p       = -1*ones(8,numel(fc_res));

for s = 1 : nSubjs
   
    pre_v  = cellfun(@get_Stats,fc_res{s}.Cpre); % pre-resection
    v1_v   = cellfun(@get_Stats,fc_res{s}.C1);   % naive virtual resection
    v2_v   = cellfun(@get_Stats,fc_res{s}.C2);   % virtual resection with partialization
    post_v = cellfun(@get_Stats,fc_res{s}.Cpost); % post-resection
    
   
    [h(1,s),p(1,s)] = kstest2(v1_v,post_v,'Tail','unequal','Alpha',alpha_level);
   
    [h(2,s),p(2,s)] = kstest2(v2_v,post_v,'Tail','unequal','Alpha',alpha_level);
  
   
    
    if(h(1,s) == 0 )
        estimation_fc_v1(s) = 1;
    end 
    if(h(2,s) == 0)
        estimation_fc_v2(s) = 1;
    end 
    
    if(any( pre_v > ref_mag+ref_vari )) % any value bigger than the reference interval
        answer_allowed(s) = 1;
         
        if(contains(fc_res{s}.so,'1A_AED'))
            [h(3,s),p(3,s)] = kstest2(pre_v,post_v,'Tail','smaller','Alpha',alpha_level);

            [h(4,s),p(4,s)] = kstest2(pre_v,v1_v,'Tail','smaller','Alpha',alpha_level);

            [h(5,s),p(5,s)] = kstest2(pre_v,v2_v,'Tail','smaller','Alpha',alpha_level);

        else
            [h(6,s),p(6,s)] = kstest2(pre_v,post_v,'Tail','larger','Alpha',alpha_level);

            [h(7,s),p(7,s)] = kstest2(pre_v,v1_v,'Tail','larger','Alpha',alpha_level);

            [h(8,s),p(8,s)] = kstest2(pre_v,v2_v,'Tail','larger','Alpha',alpha_level);
        end
        
      
    end
end



% save kstest p-value table
res_tbl = array2table(round(p,5),...
                      'VariableNames',subjNames,...
                       'RowNames',{'VRnaive_vs_post',...
                                   'VRpartialization_vs_post',...
                                   'pre vs_post_good',...
                                   'naive_vs_post_good',...
                                   'partialization_vs_post_good',...
                                   'pre vs_post_bad',...
                                   'naive_vs_post_bad',...
                                   'partialization_vs_post_bad'...
                                    ...
                                   });

writetable(res_tbl,fullfile(outFolder,tabFname),'FileType','text','Delimiter','tab','WriteRowNames',1);

fid = fopen(fullfile(outFolder,resFname),'w');

idx_good = ~idx_bad_out;

fprintf(fid,'successful estimation fc naive: %i / %i\n',sum(estimation_fc_v1),nSubjs);
fprintf(fid,'successful estimation fc partialization: %i / %i\n',sum(estimation_fc_v2),nSubjs);
fprintf(fid,'good outcome : %i / %i bad outcome  : %i / %i\n',sum(idx_good),nSubjs,sum(idx_bad_out),nSubjs);
fprintf(fid,'applicability  : %i / %i not applicability  : %i / %i\n',sum(answer_allowed),nSubjs,sum(~answer_allowed),nSubjs);
fprintf(fid,'pre > post (good): %i / %i\n',sum( h(3,answer_allowed & idx_good) == 1),sum(answer_allowed & idx_good));
fprintf(fid,'pre > naive (good): %i / %i\n',sum(h(4,answer_allowed & idx_good) == 1),sum(answer_allowed & idx_good));
fprintf(fid,'pre > partialization (good): %i / %i\n',sum(h(5,answer_allowed & idx_good) == 1),sum(answer_allowed & idx_good));
fprintf(fid,'pre < post (bad): %i / %i\n',sum(h(6,answer_allowed & idx_bad_out) == 1),sum(answer_allowed & idx_bad_out));
fprintf(fid,'pre < naive (bad): %i / %i\n',sum(h(7,answer_allowed & idx_bad_out) == 1),sum(answer_allowed & idx_bad_out));
fprintf(fid,'pre < partialization (bad): %i / %i\n',sum(h(8,answer_allowed & idx_bad_out) ==1 ),sum(answer_allowed & idx_bad_out));


fclose(fid);


if(~exist(outFolder,'dir'))
    mkdir(outFolder)
end

% all 

cfgP           = [];
cfgP.X_lim     = [0 4.5];
cfgP.Y_lim     = [0 0.5];
cfgP.nr        = 3;
cfgP.nc        = 4;
cfgP.ref_mag   = ref_mag;
cfgP.ref_vari  = ref_vari; 
cfgP.subjNames = subjNames;
cfgP.fc_res    = fc_res;

fig = plot_fig(cfgP);

% applicability 
cfgP.subjNames = subjNames(logical(answer_allowed));
cfgP.fc_res    = fc_res(logical(answer_allowed));
cfgP.nr        = 2;
cfgP.nc        = 4;
fig(2) = plot_fig(cfgP);

% not applicability
cfgP.nr        = 2;
cfgP.nc        = 4;
cfgP.subjNames = subjNames(logical(~answer_allowed));
cfgP.fc_res    = fc_res(logical(~answer_allowed));

fig(3) = plot_fig(cfgP);

if(numel(fig) == 1)
        set(fig,'WindowState','fullscreen')
        outfile = fullfile(outFolder,strcat('all_subjects','_',fc_res{1}.fc_type,'_','global'));
        print(fig,outfile,'-djpeg');
        close(fig)
else
    for s = 1 : numel(fig)

        set(fig(s),'WindowState','fullscreen')
        outfile = fullfile(outFolder,strcat('fig','_',num2str(s)));
        print(fig(s),outfile,'-djpeg');
        close(fig(s))
    end
end




function f = plot_fig(cfg)
   
   
   subjNames = cfg.subjNames;
   fc_res    = cfg.fc_res;
   X_lim     = cfg.X_lim;
   Y_lim     = cfg.Y_lim; 
   nr        = cfg.nr;
   nc        = cfg.nc; 
   ref_mag   = cfg.ref_mag;
   ref_vari  = cfg.ref_vari;
   
   f         = figure;
   nSubjs    = numel(fc_res);  
   
   
   for s = 1 : nSubjs 
    pre_v  = cellfun(@get_Stats,fc_res{s}.Cpre); % pre-resection
    v1_v   = cellfun(@get_Stats,fc_res{s}.C1);   % naive virtual resection
    v2_v   = cellfun(@get_Stats,fc_res{s}.C2);   % virtual resection with partialization
    post_v = cellfun(@get_Stats,fc_res{s}.Cpost); % post-resection

    subplot(nr,nc,s)


    val    = [pre_v , post_v, v1_v, v2_v];
    labels = [repmat({'3) pre'},1,numel(pre_v)) , repmat({'4) post'},1,numel(post_v)), ...
               repmat({'1) VR_{Naive}'},1,numel(v1_v)) , repmat({'2) VR_{partialization}'},1,numel(v2_v)), ...
             ];

    violinplot(val,labels,'showMean',true);     
    xtickangle(45) 
    hold

    ref_x  = 0:numel(unique(labels))+1;
    ref_y  = repmat(ref_mag,1,numel(unique(labels))+2);
    errBar = repmat(ref_vari,1,numel(unique(labels))+2);

    shadedErrorBar(ref_x, ref_y,errBar, 'lineprops', '-g')

    line([0 numel(unique(labels))+1],[ref_mag ref_mag],'LineStyle','--','color','k')

    ylabel('Global Connectivity')
    title(subjNames{s})

    %ylim(Y_lim)
    set(gca,'XLim',X_lim)
   end

function [avg_fc]= get_Stats(m)

        % compute average global connectivity
        get_Gfc = @(x) sum(sum(x))/(size(x,1)^2-size(x,1));
        
        % put zeros on the diagonal
        m(1:length(m)+1:end)   = 0;
        
        avg_fc = get_Gfc(m); 
        
       
       
        
        
        
        