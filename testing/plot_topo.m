% plot topoplot of the bipolar montage in two direction (horizontal and vertical) for a grid 5X4
% 
% INPUT
%
% d  - fieldtrip data structure
% V  - eigenvector matrix (channel X #eigenvectors)
 function plot_topo(d,V) 
  
    mycomp           = d;
    mycomp.topo      = V;
    mycomp.unmixing  = inv(V);
    mycomp.topolabel = d.label;
    
    nr = ceil(sqrt(length(V)));
    nc = nr;
    for j = 1 : numel(mycomp.label)
     
     subplot(nr,nc,j)
     
     cfg           = [];
     cfg.layout    = 'grid_5x4_lay_square.mat';
     cfg.component = [j];
     cfg.zlim      = 'maxabs';
     %cfg.highlight = 'labels';
     cfg.interplimits = 'electrodes';
     cfg.style     = 'straight';
   
     ft_topoplotIC(cfg, mycomp)
     colormap('jet')
     
    end