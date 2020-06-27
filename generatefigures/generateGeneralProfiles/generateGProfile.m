clear all; close all;
%% load shape: https://www.mcc.co.mercer.pa.us/dps/state_fips_code_listing.htm
points = loadGeojson('../..\ipvoronoi\data\All_states\AL_CON_Districts\PPoints\1_t4.geojson');

% normalize points
points = points - mean(points);
points = points/max(max(points));

plot(polyshape(points)); hold all; axis equal;

name = 'G';

if exist(['cached' name '.mat'],'file')
    load(['cached' name '.mat']);
else
    params.dx = .1;
    params.dc = .5;
    params.rlimlim = .1;
    params.nHeightPixels = 100;
    params.Visualize = 1;
    params.useNeckRlim = 1;
    params.debug = 0;
    params.nMorphOpenSamples = 50;
    params.postprocess = 1;
    params.postprocessclosingradius = .01;
    res = aggregateProcessing(points, params);
    save(['cached' name '.mat'],'res','params','points');
end

visparams.nticks = 7;
visparams.supressNeckShape = true;
visparams.tvscale = .115;
visparams.scale = 1.5;
visparams.tvdownshift = 3.2;
visparams.moupshift = 1;
visparams.madownshift = 1.2;
visparams.legend = false;
visparams.name = 'District 1';
fig = visualizeProfile(points, res, params, visparams);
exportgraphics(fig,[name 'Profile.pdf'],'ContentType','vector')
