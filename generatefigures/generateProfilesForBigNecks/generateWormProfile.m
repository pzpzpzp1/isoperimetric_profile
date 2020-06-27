clear all; close all;
load('wormPoints.mat');

if exist('cachedWorm.mat','file')
    load('cachedWorm.mat');
else
    params.dx = .1;
    params.dc = .1;
    params.rlimlim = .07;
    params.nHeightPixels = 100;
    params.Visualize = 1;
    params.useNeckRlim = 1;
    params.debug = 0;
    params.nMorphOpenSamples = 100;
    params.postprocess = 1;
    params.postprocessclosingradius = .01;
    res = aggregateProcessing(points, params);
    save('cachedWorm.mat','res','params','points');
end

visparams.nticks = 6;
visparams.supressNeckShape = true;
visparams.tvscale = .075;
visparams.name = 'Worm';
visparams.legend = false;
fig = visualizeProfile(points, res, params, visparams);
exportgraphics(fig,'wormProfile.pdf','ContentType','vector')

%{
regs = regions(res.MARes.popens(end));

figure; hold all; axis equal;
plot(regs(1),'facecolor','green')
plot(regs(2),'facecolor','red')
plot(regs(3),'facecolor','blue')
plot(regs(4),'facecolor','cyan','linewidth',20)

simped = rmholes(regs(1));
[area(simped ) perimeter(simped )]
[area(polyshape(points)) perimeter(polyshape(points))]

plot(polyshape(points),'facecolor','red')

%}
