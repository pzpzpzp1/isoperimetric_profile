close all; clear all;
%% load shape: https://www.mcc.co.mercer.pa.us/dps/state_fips_code_listing.htm
% points = loadGeojson('.\ipvoronoi\data\All_states\Points\10_t2_2.geojson');
% points = loadGeojson('.\ipvoronoi\data\All_states\Points\09_t3_4.geojson');
% points = loadGeojson('.\ipvoronoi\data\All_states\Points\06_t3_2.geojson');
% points = loadGeojson('.\ipvoronoi\data\All_states\Points\06_t1.geojson');
% points = loadGeojson('.\ipvoronoi\data\All_states\Points\06_t2.geojson');
% points = loadGeojson('.\ipvoronoi\data\All_states\Points\06_t4.geojson');
% points = loadGeojson('.\ipvoronoi\data\All_states\Points\05_t3_3.geojson');
% points = loadGeojson('.\ipvoronoi\data\All_states\Points\06_t5.geojson');
% points = loadGeojson('.\ipvoronoi\data\All_states\Points\01_t3_4.geojson');
% points = userinput_points('cachedshape.lines'); % generate new shape. optionally save it.
% t = linspace(0,2*pi,50)'; points = [cos(t), sin(t)]; % analytic circle
points = loadCachedShape();
% points = loadPerturbedCircle(40, .025);

%% extract some basic info
BBcenter = [min(points)+ max(points)]/2;
points = points - BBcenter;
shape = polyshape(points); 
totalArea = area(shape); 
totalPerimeter = perimeter(shape);
dc = .05;
dx = .05;
pointsPerMAEdge = 20;
nMorphOpenSamples = 100;
rlimlim = .05;
debug = 0;
Visualize = 1;
fh = visualizeMedAxis(points);

%% process morph open bound
tic; [fullperimeters0, fullareas0, keypoints, popens0, perodes0] = morphOpenBound(points, Visualize, [nMorphOpenSamples, .04]);
tmopen = toc;
rlim = max(keypoints.minRadNeck*.95, rlimlim); 
tic; [fullperimetersUpperBound2, fullareasMA2, popens2, props] = isoperim_profile_v3(points, Visualize, [dc, rlim, dx], debug);
tMA2 = toc;
startAreaPercentage = keypoints.inrArea/totalArea;

%% process TV bound with different ways of computing discrete TV
nHeightPixels = 100;
MidAreaSampleNumber = 13;
EndAreaSampleNumber = 7;
percentAreasMid = linspace(startAreaPercentage, keypoints.areaMinNeck/keypoints.areaMax, MidAreaSampleNumber);
percentAreasEnd = linspace(keypoints.areaMinNeck/keypoints.areaMax, .9999, EndAreaSampleNumber);
percentAreas = [percentAreasMid percentAreasEnd(2:end)];
tic; [fullperimetersLowerBound1, fullareasTV1, images1] = getTVProfile(points, nHeightPixels, 'cvx', percentAreas);
tTVcvx = toc;
tic; [fullperimetersLowerBound2, fullareasTV2, images2] = getTVProfile(points, nHeightPixels, 'admm', percentAreas);
tTVadmm = toc;

%% times
times = [tmopen, tMA2, tTVcvx, tTVadmm];

%% Visualize results
figure; hold all;
plot(fullareas0, fullperimeters0, 'g-','linewidth',2)
plot(fullareasMA2, fullperimetersUpperBound2, 'b-','linewidth',2)
plot(fullareasTV1, fullperimetersLowerBound1, 'm-','linewidth',1)
plot(fullareasTV2, fullperimetersLowerBound2, 'c-','linewidth',1)

xline(keypoints.inrArea,'k','linewidth',2)
xline(keypoints.areaMinNeck,'g')
xline(keypoints.areaMax,'k','linewidth',2)
for i=1:size(keypoints.eventAreas,1)
    xline(keypoints.eventAreas(i,1),'b')
    xline(keypoints.eventAreas(i,2),'r')
end
legend('morph open','MA1','MA2','TV-3','TV');

