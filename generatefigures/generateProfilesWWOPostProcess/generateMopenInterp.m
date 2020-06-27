clear all; close all;
load('bigbarbell.mat');

if exist('cachedBB0.mat','file')
    load('cachedBB0.mat');
else
    params.dx = .01;
    params.dc = 1000;
    params.rlimlim = .07;
    params.nHeightPixels = 5;
    params.Visualize = 1;
    params.useNeckRlim = 1;
    params.debug = 0;
    params.nMorphOpenSamples = 200;
    params.postprocess = 1;
    params.postprocessclosingradius = .01;
    params.skipTV = 1;
    res = aggregateProcessing(points, params);
    save('cachedBB0.mat','res','params','points');
end

fh=figure; ax0 = gca; hold all; set(gcf,'color','w');
fh.Renderer = 'Painters';
ylabel('Perimeter'); xlabel('Area'); title('Morph Open Interpolation');
fh.OuterPosition = [1   1    380  380];
plot(res.ORes.areas, res.ORes.perims, 'b-','linewidth',2)
xlim([0.8 1.85]); 
ylim([2.5 6]);
drawnow;


area1 = res.ORes.keypoints.eventAreas(1,1);
area2 = res.ORes.keypoints.eventAreas(1,2);
xline(area1,'g','linewidth',2);
xline(area2,'r','linewidth',2);
erodedview = morpherode(polyshape(points),.52);
preerodedview = morpherode(polyshape(points),.525);
drawnow;

mockBBc = (min(points) + max(points))/2;
shape = polyshape(points);
areas = linspace(area1,area2,4);
anisoscaleORat = fh.OuterPosition;
anisoscaleO = anisoscaleORat(4)/anisoscaleORat(3);
anisoscaleRat = get(gca,'DataAspectRatio');
anisoscale = anisoscaleRat(2)/anisoscaleRat(1)/anisoscaleO;    
scale = .38;
for i=1:numel(areas)
    areaI = areas(i);
    [~, ind] = min(abs(res.ORes.areas-areaI));
    center = [areaI 3];
    shapei = res.ORes.popens(ind);
    shapei.Vertices = shapei.Vertices;
    
    oshape = transformShape(shape, scale, center, anisoscale);
    plot(oshape,'facecolor','blue','linewidth',1)

    transshape = transformShape(shapei, scale, center, anisoscale, mockBBc);
    plot(transshape,'facecolor',[.5 1 .5],'facealpha',1)

    if i==1
        erodeshape = preerodedview;
    else
        erodeshape = erodedview;
    end
    erodeshape = transformShape(erodeshape, scale, center, anisoscale, mockBBc);
    plot(erodeshape,'facecolor',[1 .5 .5],'facealpha',1,'edgecolor','r')
end
exportgraphics(fh,'mopenInterp.pdf','ContentType','vector')


