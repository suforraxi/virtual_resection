function pre_vs_post_main()
addpath('/home/matteo/Desktop/git_rep/notBoxPlot/code/')
addpath('/home/matteo/Desktop/git_rep/Violinplot-Matlab/')
load('/home/matteo/Desktop/virtual_resection/tables/summary_tbl_coh.mat')
%close all 

 %rTbl = rTbl(~(rTbl.biomarker>0.6),:)
 %figure
 %plot_viol(rTbl)
 
 %remove outliers
 idx        = true(size(rTbl,1),1);
 idx([8,22])= false;
 rTbl       = rTbl(idx,:);
 
 uSubjName = unique(rTbl.subjName);
 nSubj     = numel(uSubjName);
 nr        = ceil(sqrt(nSubj));
 
 figure
 
 for i = 1 : nSubj
    
    subplot(nr,nr,i)
    
    subj_T = rTbl(strcmp(rTbl.subjName,uSubjName{i}),:);
    if(any(strcmp(subj_T.sitType,'Pre')) && any(strcmp(subj_T.sitType,'Post')) )
        plot_subj(subj_T)
    end
    
 end
 
 figure
 post_T = rTbl(strcmp(rTbl.sitType,'Post'),:);
 
 good_idx = ~cellfun(@isempty,regexp(post_T.seizOut,'1(A|B)\w*'));
 bad_idx  = ~good_idx;
 a        = post_T.biomarker(good_idx);
 b        = post_T.biomarker(bad_idx);

 plot_a_vs_b(a,b);
 
 
 

 
 function plot_subj(rTbl)
 
 post_val = rTbl.biomarker(strcmp('Post',rTbl.sitType));  
 pre_val  = rTbl.biomarker(strcmp('Pre',rTbl.sitType));
 virt_val = rTbl.v_biomarker(strcmp('Pre',rTbl.sitType));
 %[h,p]    = kstest2(post_val,pre_val,'Tail','larger');

% catNames = [repmat({'PRE'},size(pre_val),1 ); repmat({'POST'},size(post_val),1)] ;
%   
%   viol = violinplot([pre_val;post_val], catNames);
%   viol(1).ViolinColor = [0 0 1];
%   viol(2).ViolinColor = [1 0 0];

catNames = [repmat({'PRE'},size(pre_val),1 ); ...
            repmat({'POST'},size(post_val),1); ...
            repmat({'VIRT'},size(virt_val),1);
            ] ;
   
%viol                = violinplot([pre_val;post_val;virt_val], catNames);

bar([pre_val;post_val;virt_val]);
xticklabels(catNames)
 %g = [repmat(0,size(pre_val),1 ); repmat(1,size(post_val),1)] ;
 %notBoxPlot([pre_val;post_val],g)
 ylim([0 1])
 ylabel('Coherence')
 
 
 
 if(numel(unique(rTbl.subjName))==1)
     sName = rTbl.subjName{1};
     so    = rTbl.seizOut{1};
     title(sprintf('%s %s',sName,so),'Interpreter','none')
 end    

 function plot_a_vs_b(a,b)
     
   
 [h,p]    = kstest2(a,b,'Tail','unequal');

% catNames = [repmat({'PRE'},size(pre_val),1 ); repmat({'POST'},size(post_val),1)] ;
%   
%   viol = violinplot([pre_val;post_val], catNames);
%   viol(1).ViolinColor = [0 0 1];
%   viol(2).ViolinColor = [1 0 0];

 g = [repmat(0,size(a),1 ); repmat(1,size(b),1)] ;
 notBoxPlot([a;b],g)
 ylim([0 1])
    

 title(sprintf('p: %.4f',p),'Interpreter','none')
 
 
function plot_pre_post_virt(a,b,c)


%viol(1).ViolinColor = [0 0 1];
%viol(2).ViolinColor = [1 0 0];

ylim([0 1])
     
         
 
     
         