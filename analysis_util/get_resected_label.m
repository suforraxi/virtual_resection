
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