clear all; close all;
load('bigbarbell.mat');

if exist('cachedBB.mat','file')
    load('cachedBB.mat');
else
    params.dx = .1;
    params.dc = .1;
    params.rlimlim = .07;
    params.nHeightPixels = 5;
    params.Visualize = 1;
    params.useNeckRlim = 1;
    params.debug = 0;
    params.nMorphOpenSamples = 100;
    params.postprocess = 1;
    params.postprocessclosingradius = .01;
    params.skipTV = 1;
    res = aggregateProcessing(points, params);
    save('cachedBB.mat','res','params','points');
end

if exist('cachedBB2.mat','file')
    load('cachedBB2.mat');
else
    params2.dx = .1;
    params2.dc = .1;
    params2.rlimlim = .07;
    params2.nHeightPixels = 5;
    params2.Visualize = 1;
    params2.useNeckRlim = 1;
    params2.debug = 0;
    params2.nMorphOpenSamples = 100;
    params2.postprocess = 0;
    params2.postprocessclosingradius = .01;
    params2.skipTV = 1;
    res2 = aggregateProcessing(points, params2);
    save('cachedBB2.mat','res2','params2','points');
end

%% clean up data. Exclude TV nans and exclude past rlimlim
% show MA data up to rlim
rlimused = res.ORes.keypoints.minRadNeck;
rlimArea = area(morphopen(polyshape(points), rlimused));
keepindsMA = res.MARes.areas < rlimArea;

fh=figure; ax0 = gca; hold all; set(gcf,'color','w');
fh.Renderer = 'Painters';
ylabel('Perimeter'); xlabel('Area'); title('IP Postprocessing');
fh.OuterPosition = [1   1    355.2000  408.0000];
plot(res.MARes.areas(keepindsMA), res.MARes.perims(keepindsMA), 'r-','linewidth',2)
plot(res2.MARes.areas(keepindsMA), res2.MARes.perims(keepindsMA), 'g-','linewidth',2)
xlim([res.ORes.keypoints.inrArea*.95 rlimArea*1.05]); 
xline(res.ORes.keypoints.inrArea,'k','linewidth',2,'label','Max Inscribed Circle')
xline(rlimArea,'k','linewidth',2,'label','Minimal Neck','labelverticalalignment','bottom','labelhorizontalalignment','left')
ylim([2*pi*res.ORes.keypoints.inr *.95, max(res.MARes.perims(keepindsMA))*1.05]);
legend('Postprocessed','Raw')
exportgraphics(fh,'PostprocessingvsRaw.pdf','ContentType','vector')

%% generate data send to justin
shapePre = res2.MARes.popens(end-7);
shapePost = postProcessCandidatePoly(polyshape(points), shapePre, .01);
save('pre_and_post.mat','shapePre','shapePost');



