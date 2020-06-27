clear all; close all;
load('bigbarbell.mat');

if exist('cachedBB.mat','file')
    load('cachedBB.mat');
else
    params.dx = .1;
    params.dc = .01;
    params.rlimlim = .07;
    params.nHeightPixels = 100;
    params.Visualize = 1;
    params.useNeckRlim = 1;
    params.debug = 0;
    params.nMorphOpenSamples = 100;
    params.postprocess = 1;
    params.postprocessclosingradius = .01;
    res = aggregateProcessing(points, params);
    save('cachedBB.mat','res','params','points');
end

visparams.nticks = 5;
visparams.tvscale = .075;
visparams.name = 'Hourglass'
visparams.legend = 0;
fig=visualizeProfile(points, res, params, visparams);
exportgraphics(fig,'bbProfile.pdf','ContentType','vector')


