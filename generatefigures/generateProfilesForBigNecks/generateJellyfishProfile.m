clear all; close all;
load('jellyfish.mat');

if exist('cachedJellyfish.mat','file')
    load('cachedJellyfish.mat');
else
    params.dx = .3;
    params.dc = .5;
    params.rlimlim = .07;
    params.nHeightPixels = 100;
    params.Visualize = 1;
    params.useNeckRlim = 1;
    params.debug = 0;
    params.nMorphOpenSamples = 100;
    params.postprocess = 1;
    params.postprocessclosingradius = .01;
    res = aggregateProcessing(points, params);
    save('cachedJellyfish.mat','res','params','points');
end

visparams.nticks = 10;
visparams.tvscale = .075;
visparams.name = 'Jellyfish';
visparams.legend = true;
fig=visualizeProfile(points, res, params, visparams);
exportgraphics(fig,'jellyfishProfile.pdf','ContentType','vector')

%% debugging stuff
%{
bs = polyshape(points);
nareas = [];
nperims = [];
for i =1:numel(res.MARes.popens)
    poly = res.MARes.popens(i);
    poly = postProcessCandidatePoly(bs, poly, .01)
    nareas(i) = area(poly);
    nperims(i) = perimeter(poly);
end

figure; hold all;
plot(res.MARes.areas, res.MARes.perims, 'm.-','linewidth',1)
plot(nareas, nperims,'c')
yline(perimeter(bs))
%}
