clear; close all;

% points = userinput_points(); save('pointsopen.mat','points');
% points = loadPerturbedCircle(3, .025);
load starPoints.mat;

%% show dilate and erode
rad = .6;
ps = polyshape(points);

%% morphopen
%rads = [.1 .3 .5 .7];
rads = [.175 .35 .525 .7]/2;
popen = polyshape.empty(0,numel(rads));
for i=1:numel(rads)
    popen(i) = morphopen(ps, rads(i));
end

f2 = figure; axis equal; hold all; set(gcf,'color','w'); axis off;
plot(ps,'facecolor','blue','facealpha',.06,'linewidth',2);
for i=1:numel(rads)
    plot(subtract(ps, popen(i)),'facecolor','blue','facealpha',.06);
end
exportgraphics(f2, 'morphStarOpen.pdf', 'contenttype','vector')
