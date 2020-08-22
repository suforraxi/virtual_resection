% 
% Setting path for all the dependencies for virtual resection project
% All the matlab functions and toolboxes required
% They are freely available on the web


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


% git root
git_root = '/Users/matte/Desktop/git_rep/';

% functions to import BIDS structures for ioECoG  (ithub.com/suforraxi/ieeg_respect_bids)         
addpath(fullfile(git_root,'ieeg_respect_bids'))
addpath(fullfile(git_root,'ieeg_respect_bids','external'))
addpath(fullfile(git_root,'ieeg_respect_bids','importBIDS'))
addpath(fullfile(git_root,'ieeg_respect_bids','trc2bids'))
addpath(fullfile(git_root,'ieeg_respect_bids','micromed_utils'))


% fieldtrip (https://github.com/fieldtrip)
addpath((fullfile(git_root,'fieldtrip')));
addpath((fullfile(git_root,'fieldtrip','fieldtrip_private')));


% json functions (https://github.com/fangq/jsonlab)
addpath(fullfile(git_root,'jsonlab')) 


% virtual resection
addpath(fullfile(git_root,'virtual_resection'))
%addpath(fullfile(git_root,'virtual_resection','import_trc'))
addpath(fullfile(git_root,'virtual_resection','montage'))
addpath(fullfile(git_root,'virtual_resection','analysis_util'))
addpath(fullfile(git_root,'virtual_resection','external'))

% add violin plot functions (https://github.com/bastibe/Violinplot-Matlab)
addpath('/Users/matte/Desktop/git_rep/Violinplot-Matlab/')

% add shaded error bar functions (https://github.com/raacampbell/shadedErrorBar)
addpath('/Users/matte/Desktop/git_rep/shadedErrorBar/')


