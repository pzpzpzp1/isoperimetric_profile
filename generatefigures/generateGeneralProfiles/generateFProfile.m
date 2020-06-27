clear all; close all;
load mozpoints.mat;
points = points - mean(points);
points = points/max(max(points));
[x y] = polyshape(points).boundary;
points = [x y];
name = 'F';

if exist(['cached' name '.mat'],'file')
    load(['cached' name '.mat']);
else
    params.dx = .02;
    params.dc = 3.5;
    params.rlimlim = .001;
    params.nHeightPixels = 100;
    params.Visualize = 1;
    params.useNeckRlim = 1;
    params.debug = 0;
    params.nMorphOpenSamples = 200;
    params.postprocess = 1;
    params.postprocessclosingradius = .01;
    res = aggregateProcessing(points, params);
    save(['cached' name '.mat'],'res','params','points');
end

visparams.nticks = 7;
visparams.supressNeckShape = true;
visparams.tvscale = .13;
visparams.scale = 1.543/2;
visparams.tvdownshift = 2.4;
visparams.moupshift = 1;
visparams.madownshift = 1;
visparams.legend = false;
visparams.maxinrdir = 'left';
visparams.name = 'Mozambique';
fig = visualizeProfile(points, res, params, visparams);
exportgraphics(fig,[name 'Profile.pdf'],'ContentType','vector')
