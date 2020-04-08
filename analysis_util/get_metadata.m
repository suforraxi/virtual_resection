% function to obtain metadata regarding the recordings
% if a channel was resected or not during surgery
% if a channel contains any artefact during the recordings
% 
% INPUT
% dataDir  - root folder for raw data (data should be organized in BIDS)
% fileName - filename for which to read the metadata (.vhdr) 
%
% OUTPUT
% res_channel - cell array with the name of the resected channels
%
% artefact_T  - table with the following variables
%               
%               type    name of the artefact (we have just generic artefact called 'artefact')        
%               start   sample time of the original data when the artefact starts        
%               stop    sample time of the original data when the artefact ends
%               channel channel name where the artefact was found
 
function [ res_channel, artefact_T] = get_metadata(dataDir,fileName)


[~,fileName,ext] = fileparts(fileName);
subFolder        = split(fileName,'_'); 

annotFile        = fullfile(dataDir,...
                     subFolder{1},...
                     subFolder{2},...
                     'ieeg',...
                      strcat(fileName,'_annotations.tsv'));


tsv_annots       = readtable(annotFile, 'Delimiter', 'tab', 'FileType', 'text', 'ReadVariableNames', true);
idx_resected     = strcmp(tsv_annots.type,'resected');
idx_edge         = strcmp(tsv_annots.type,'edge');

if(~isempty(idx_resected))
    res_channel      = tsv_annots.channel(idx_resected,1);
else
    res_channel      = [];
end
%edge_channel     = tsv_annots.channel(idx_edge,1);

idx_artefacts    = strcmp(tsv_annots.type,'artefact');
if(~isempty(idx_artefacts))
    artefact_T        = tsv_annots(idx_artefacts,:);
else
    artefact_T        = [];
end
