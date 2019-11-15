%% remove linear contribution from resetcted channels

function [nResData,resData] = remove_resected_influence(cfg,data)



[ res_channel, artefact_T] = get_metadata(cfg.bidsFolder,cfg.sitFname);
 res_notR_cut              = cell(size(data.label));
 
 for ch = 1 : numel(data.label)
    res_notR_cut{ch,1} = get_resected_label(data.label{ch},res_channel);
 end
 
 
 
nResData = data;
resData  = data;

nTrials = numel(data.trial);

idx_res = cellfun(@isempty,regexp(res_notR_cut,'NRES'));

for t = 1 : nTrials

    m = data.trial{t};

    X = m(idx_res,:);
    Y = m(~idx_res,:);

    nchTarget = size(Y,1); 
    reg       = X';
    R         = [];
    for ch = 1 : nchTarget
        target  = Y(ch,:)';
        [B,~,R(ch,:)] = regress(target,reg);
    end


    nResData.label    = data.label(~idx_res);
    nResData.trial{t} = R;
    
    resData.label    = data.label(idx_res);
    resData.trial{t} = data.trial{t}(idx_res,:);
    
end