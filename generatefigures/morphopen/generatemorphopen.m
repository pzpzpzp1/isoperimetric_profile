clear; close all;

% points = userinput_points(); save('pointsopen.mat','points');
% points = loadPerturbedCircle(3, .025);
load pointsopen.mat;
R = [0 -1; 1 0;]; points = points*R;

%% show dilate and erode
rad = .6;
ps = polyshape(points);
[popen, perode] = morphopen(ps, rad);
[pclose, pdilate] = morphclose(ps, rad);

f1 = figure; axis equal; hold all; set(gcf,'color','w'); axis off;
plot(pdilate,'facecolor',[.2 .8 .2],'facealpha',1,'edgecolor','none');
plot(ps,'facecolor',[.7 .7 1],'facealpha',1,'linewidth',2);
plot(perode,'facecolor',[.8 .3 .3],'facealpha',1,'linewidth',1,'edgecolor','none');
exportgraphics(f1, 'morphA.pdf', 'contenttype','vector')

%% morphopen
%rads = [.1 .3 .5 .7];
rads = [.175 .35 .525 .7]*2;
popen = polyshape.empty(0,numel(rads));
for i=1:numel(rads)
    popen(i) = morphopen(ps, rads(i));
end

f2 = figure; axis equal; hold all; set(gcf,'color','w'); axis off;
plot(ps,'facecolor','blue','facealpha',.06,'linewidth',2);
for i=1:numel(rads)
    plot(subtract(ps, popen(i)),'facecolor','blue','facealpha',.06);
end
exportgraphics(f2, 'morphB.pdf', 'contenttype','vector')

%% morphclose
rads = [.4 .8 1.2 1.6 2];
pclose = polyshape.empty(0,numel(rads));
for i=1:numel(rads)
    pclose(i) = morphclose(ps, rads(i));
end

f3 = figure; axis equal; hold all; set(gcf,'color','w'); axis off;
plot(ps,'facecolor','blue','facealpha',.06,'linewidth',2);
for i=1:numel(rads)
    plot(pclose(i),'facecolor','blue','facealpha',.06);
end
exportgraphics(f3, 'morphC.pdf', 'contenttype','vector')