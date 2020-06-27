clear all; close all;
%% load shape: https://www.mcc.co.mercer.pa.us/dps/state_fips_code_listing.htm

points = loadGeojson('../..\ipvoronoi\data\All_states\Points\06_t3.geojson');

% normalize points
points = points - mean(points);
points = points/max(max(points));

plot(polyshape(points)); hold all; axis equal;

name = 'B';

if exist(['cached' name '.mat'],'file')
    load(['cached' name '.mat']);
else
    params.dx = .1;
    params.dc = 1;
    params.rlimlim = .07;
    params.nHeightPixels = 100;
    params.Visualize = 1;
    params.useNeckRlim = 1;
    params.debug = 0;
    params.nMorphOpenSamples = 100;
    params.postprocess = 1;
    params.postprocessclosingradius = .01;
    res = aggregateProcessing(points, params);
    save(['cached' name '.mat'],'res','params','points');
end

visparams.nticks = 10;
visparams.supressNeckShape = true;
visparams.tvscale = .1;
%visparams.scale = .7143/2;
%visparams.tvdownshift = 1.7;
%visparams.moupshift = .6;
%visparams.madownshift = .7;
visparams.name = 'California';
visparams.legend = false;
fig = visualizeProfile(points, res, params, visparams);
exportgraphics(fig,[name 'Profile.pdf'],'ContentType','vector')
