% Get resected channels for bipolar montage
% 
% INPUT
%      chName      - cell with the label for bipolar derivation (channels)
%      res_channel - cell with the labels with resected channels as
%                    recorded (common reference)
% OUTPUT
%
%     res_notR_cut - string indicating if the bipolar channel was
%                    resected (RES) both monopolar channels resected 
%                    not resected (NRES) both monopolar channels not resected  
%                    cut (CUT) one monopolar channel resected the other not
%                    resected
%
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

function res_notR_cut = get_resected_label(chName,res_channel)



ch = regexpi(chName,'(\w*\d+)[N]?-[N]?(\w*\d+)','tokens');

ch1 = ch{1}{1}; 
ch2 = ch{1}{end};


res_notR_cut = 'NRES';
if( any(strcmp(ch1,res_channel)) || any(strcmp(ch2,res_channel)) )
    if(any(strcmp(ch1,res_channel)) && any(strcmp(ch2,res_channel)))
        res_notR_cut = 'RES';
    else
        res_notR_cut = 'CUT';
    end
end