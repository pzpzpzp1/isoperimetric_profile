clear all; close all;
points = loadCachedShape('./candy.lines');
BBcenter = [min(points)+ max(points)]/2;
points = points - BBcenter;


if exist('cachedCandy.mat','file')
    load('cachedCandy.mat');
else
    params.dx = .1;
    params.dc = .1;
    params.rlimlim = .14;
    params.nHeightPixels = 200;
    params.Visualize = 1;
    params.useNeckRlim = 0;
    params.debug = 0;
    params.nMorphOpenSamples = 100;
    params.postprocess = 1;
    params.postprocessclosingradius = .01;
    params.skipTV = 1;
    res = aggregateProcessing(points, params);
    save('cachedCandy.mat','res','params','points');
end

% mafig = visualizeMedAxis(points);
% exportgraphics(mafig,'candyMedAxis.pdf','ContentType','vector');

visparams.nticks = 10;
visparams.tvscale = .06;
visparams.supressNeckShape = true;
visparams.name = 'Candy';
visparams.OuterPosition = [1 1 1000 400];
visparams.lgdpos = [0.3555 0.7490 0.1721 0.1627];
visparams.skiptv = 1;
%visparams.scale = .7143/2;
%visparams.tvdownshift = 1.7;
%visparams.moupshift = .6;
%visparams.madownshift = .7;

shape = polyshape(points);
rlimused = params.rlimlim;
rlimArea = area(morphopen(shape, rlimused));
keepindsMA = res.MARes.areas < rlimArea;
rnArea = res.ORes.keypoints.areaMinNeck;
rt = max(res.MARes.props.maxneckrad, res.MARes.props.maxCirc2rad/2);
inxy = res.ORes.keypoints.inxy;
leftsidetightT = getArea1(shape, rt, inxy);
fh=figure; ax0 = gca; hold all; set(gcf,'color','w');
fh.Renderer = 'Painters';
ylabel('Perimeter'); xlabel('Area'); 
title(visparams.name);
fh.OuterPosition = [1   1    700  350.2000];
bump = sqrt(var(res.ORes.perims))/2;
xlim([res.ORes.keypoints.inrArea*.9 res.ORes.keypoints.areaMax*1.05]); 
ylim([4   32.1910]);
plot(res.MARes.areas(keepindsMA), res.MARes.perims(keepindsMA), 'b-','linewidth',2)
xline(res.ORes.keypoints.inrArea,'k','linewidth',2,'label','Max inscribed circle','labelhorizontalalignment','right')
xline(res.ORes.keypoints.areaMax,'k','linewidth',2,'label','Max area','labelverticalalignment','bottom')
xline(rnArea,'k','linewidth',2,'label','Minimal Neck','labelverticalalignment','bottom','labelhorizontalalignment','left')
xline(rnArea,'m','linewidth',2)
xline(leftsidetightT,'k','linewidth',2,'label','Conservatively tight','labelhorizontalalignment','right')
xline(leftsidetightT,'m','linewidth',2)
maxarea = area(polyshape(points));

ptc = patch('vertices',[0,0,-20;leftsidetightT 0,-20; leftsidetightT 1000,-20; 0 1000,-20],'faces',[1 2 3 4],'facecolor','green'); alpha(ptc,.4);
ptc = patch('vertices',[rnArea,0,-20;...
maxarea 0,-20;...
 maxarea 1000,-20;...
 rnArea 1000,-20],'faces',[1 2 3 4],'facecolor','green'); alpha(ptc,.4);

significantAreas = [linspace(res.ORes.keypoints.inrArea, rlimArea, visparams.nticks)];
scale = .5;
anisoscaleORat = fh.OuterPosition;
anisoscaleO = anisoscaleORat(4)/anisoscaleORat(3);
anisoscaleRat = get(gca,'DataAspectRatio');
anisoscale = (anisoscaleRat(2)/anisoscaleRat(1)) / anisoscaleO;
xlimb = xlim; ylimb = ylim;
iterlist = 1:numel(significantAreas)-1;
tvdownshift = 2;
moupshift = 1;
madownshift = 1;
for i=iterlist
    areaI = significantAreas(i);

    [~,ind1] = min(abs(res.MARes.areas - areaI));
    ps1 = res.MARes.popens(ind1);
    perim1 = res.MARes.perims(ind1);
    ps1T = transformShape(ps1, scale, [areaI, perim1-madownshift*bump], anisoscale);
    plot(ps1T,'facecolor','blue','facealpha',.8)

end
ps1T = transformShape(res.MARes.popens(end), scale, [res.ORes.keypoints.areaMax, res.MARes.perims(end-1)-madownshift*bump], anisoscale);
plot(ps1T,'facecolor','blue','facealpha',.8)
xlim(xlimb); ylim(ylimb);
set(gca, 'Color', 'None')
lgd = legend('Med Axis')
lgd.Position = [0.4043    0.7894    0.2638    0.0873];
fh.OuterPosition = [267.4000  350.6000  648.8000  387.2000];
exportgraphics(fh,'candyProfile_minimal.pdf','ContentType','vector');

