% input data folder
cfg_main.bidsFolder = '/Users/matte/Desktop/git_rep/bids_4_sharing_acute/output/virtual_resection_data/data/';
% output folder where to save results
%cfg_main.outFolder  = './output/';

% subject information table (containing seizure outcome and other info)
cfg_main.info_F = './info/info_virtual_resection.tsv';


cfg_main.fc_type         = 'h2';  % functional connectivity measure
cfg_main.montage         = 'avg'; % montage to apply
    
cfg_main.lastMinute      = 0;     % flag to control the if using the last minute or the whole recordings
cfg_main.L               = 25;    % last minute
cfg_main.overlap         = 0;     % no overlapping 

cfg_main.Tlength         = 10;     % seconds of new trials
cfg_main.Toverlap        = 0;     % no overlapping between the epochs

cfg_main.notch           = [50 100 150]; % notch filtering 
%cfg_main.bpfreq          = [5  150];     % frequency interval of interest 

cfg_main.notchBool       = 'yes';
cfg_main.bpBool          = 'yes';

cfg_main.EnvType         = 'filt_hilbert_env'; % compute the envelope using hilbert
cfg_main.FConEnv         = 'yes';              % flag to control if using envelopes of the signal or raw signal

band_oi    = {[5 150],[5 15],[15 25],[30 40],[95 105],[105 150]...
             };
          
band_names = {'broadband_10s','alpha_theta','beta','low_gamma','high_gamma','very_high'...
             };

for b = 1 : 1%numel(band_oi)
    
    cfg_main.bpfreq     = band_oi{b}; 
    cfg_main.outFolder  = fullfile('.','output','bands_analysis',band_names{b});
    
    virtual_resection_main(cfg_main);
    
end