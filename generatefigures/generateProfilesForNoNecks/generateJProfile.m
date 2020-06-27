clear all; close all;
load('J.mat');

if exist('cachedJ.mat','file')
    load('cachedJ.mat');
else
    params.dx = .1;
    params.dc = .1;
    params.rlimlim = .11;
    params.nHeightPixels = 100;
    params.Visualize = 1;
    params.useNeckRlim = 0;
    params.debug = 0;
    params.nMorphOpenSamples = 100;
    res = aggregateProcessing(points, params);
    save('cachedJ.mat','res','params','points');
end
params.useNeckRlim = false;
visparams.nticks = 10;
visparams.tvscale = .06;
visparams.name = 'no-neck-J';
visparams.legend = false;
visparams.supressNeckShape = 1;
fig = visualizeProfile(points, res, params, visparams);
exportgraphics(fig,'JProfile.pdf','ContentType','vector');
