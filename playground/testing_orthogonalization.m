% testing orthogonalization



t  = 1 : 1000;
fr = [11 16 25 ];

x = zeros(numel(fr),numel(t));

for i = 1 : numel(fr)
    x(i,:) = sin(2*pi*fr(i)*t./length(t));
end

y = zeros(numel(fr),numel(t));
for i = 1 : 2
    y(i,:) = x(i,:) + x(end,:);
end
y(end,:) = x(end,:);

% orthogonalization
om = get_ortho_matrix(y,3);

% plot signals in different lines and correlation matrices
% signals
f(1) = figure;
ay_lim = [-4 25];
subplot(1,3,1)
xx = x + repmat([0, 10, 20]',1,size(x,2));
plot(xx')
title('Original')
ylim(ay_lim)
subplot(1,3,2)
yy = y + repmat([0, 10, 20]',1,size(y,2));
plot(yy')
title('Mixed')
ylim(ay_lim)

omm = om + repmat([0, 10]',1,size(om,2));
subplot(1,3,3)
plot(omm')
title('Orthogonalized')
ylim(ay_lim)

% correlation matrices
figure;
imLim = [-1 1];

subplot(2,2,1)
imagesc(corr(x'),imLim)
title('Original')

subplot(2,2,2)
mixed = corr(y');
imagesc(corr(y'),imLim)
title('Mixed')

subplot(2,2,3)
imagesc(mixed(1:2,1:2),imLim)
title('Naive')

subplot(2,2,4)

imagesc(corr(om'),imLim)
title('Orthogonalized')

colorbar('manual','Position',[0.93,0.1,0.03,0.8])



