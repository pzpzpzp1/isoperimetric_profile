clear all; close all;
%% load shape: https://www.mcc.co.mercer.pa.us/dps/state_fips_code_listing.htm


% normalize points
load Epoints.mat;
points = points - mean(points);
points = points/max(max(points));

% plot(polyshape(points)); hold all; axis equal;

name = 'E';

if exist(['cached' name '.mat'],'file')
    load(['cached' name '.mat']);
else
    params.dx = .1;
    params.dc = .16;
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

visparams.nticks = 5;
visparams.supressNeckShape = true;
visparams.tvscale = .08;
visparams.scale = 1.8143/2;
visparams.tvdownshift = 2.1;
visparams.moupshift = 1;
visparams.madownshift = .7;
visparams.legend = false;
visparams.name = 'Boxes';
fig = visualizeProfile(points, res, params, visparams);
exportgraphics(fig,[name 'Profile.pdf'],'ContentType','vector')
