
% input folder with the trc
trc_input_dir = '/home/matteo/Desktop/virtual_resection/sel_TRC/'; 
trcFiles      = dir(fullfile(trc_input_dir,'*.TRC'));

cfg          = [];
cfg.proj_dir = '/home/matteo/Desktop/new_BIDS_import/virtual_resection/';  % folder to store bids files


for i = 1 : numel(trcFiles)

    cfg.filename = fullfile(trcFiles(i).folder,trcFiles(i).name);
    
    
    [status,msg,output] = annotatedTRC2bids(cfg);
      
        

end