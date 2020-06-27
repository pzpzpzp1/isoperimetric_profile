clear all; close all;
load('starPoints.mat');

if exist('cachedStar.mat','file')
    load('cachedStar.mat');
else
    params.dx = .1
    params.dc = .1
    params.rlimlim = .11;
    params.nHeightPixels = 100;
    params.Visualize = 0;
    params.useNeckRlim = 0;
    params.debug = 0;
    params.nMorphOpenSamples = 100;
    res = aggregateProcessing(points, params);
    save('cachedStar.mat','res','params','points');
end
visparams.nticks = 10;
visparams.tvscale = .1;
visparams.supressNeckShape = 0;
visparams.name = 'Star';
visparams.legend = true;
fig = visualizeProfile(points, res, params, visparams);
exportgraphics(fig,'starProfile.pdf','ContentType','vector');

