
% input folder with the trc
trc_input_dir = '/home/matteo/Desktop/virtual_resection/sel_TRC/'; 
trcFiles      = dir(fullfile(trc_input_dir,'*.TRC'));

cfg          = [];
cfg.proj_dir = '/home/matteo/Desktop/virtual_resection/BIDS_data/';            % folder to store bids files


for i = 1 : numel(trcFiles)

    cfg.filename = fullfile(trcFiles(i).folder,trcFiles(i).name);
    
    
    [status,msg,output] = annotatedTRC2bids(cfg);
      
        

end