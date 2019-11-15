%% create a layout for a 5 by 4 grid

w = 1.4;
h = 0.7;
[x1,y1] = ndgrid(2:2:9,1:2:7);

[x2,y2] = ndgrid(1:2:9,2:2:7);

x2 = x2';
y2 = y2';

x = [ x1(:) ; x2(:)];
y = [ y1(:) ; y2(:)];

layout.pos   = [x y];
layout.label = m_data.label;%split(int2str(1:31)) ;

layout.width   = repmat(w,size(m_data.label,1),1);
layout.height  = repmat(h,size(m_data.label,1),1);

save('grid_5x4_lay.mat','layout')
%%
cfg         = [];
cfg.layout  = 'grid_5x4_lay.mat';
cfg.output  = 'grid_5x4_lay_square.mat';
cfg.outline = 'square';
cfg.mask    = 'square';

layout = ft_prepare_layout(cfg)