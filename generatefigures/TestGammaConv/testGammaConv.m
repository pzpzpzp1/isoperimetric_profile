clear all; close all;
circ = circlepoly_v3(1,[0,0],1000);
idealline = [0 0; area(circ) perimeter(circ)];


% pixels = [10:10:50, 100:100:500];
pixels = [10:10:50, 100:100:200];
percentAreas = [.02 .1:.3:.9 .98];
for i=1:numel(pixels)
    npixels = pixels(i);
	[perimeters{i}, areas{i}, images{i}] = getTVProfile(circ.Vertices, npixels, 'cvx', percentAreas);
    [perimeters_admm{i}, areas_admm{i}, images_admm{i}] = getTVProfile(circ.Vertices, npixels, 'admm', percentAreas);
end

fh=figure; %fh.WindowStyle='docked';
subplot(1,2,1); hold all;
title('DTV profile CVX fixed under mesh refinement');
plot(idealline(:,1),idealline(:,2),'linewidth',2,'color','g');
for i=1:numel(pixels)
    plot(areas{i}, perimeters{i},'.-','color', [i 0 0]/numel(pixels),'linewidth',1);
end
legend(['correct';split(num2str(pixels))]);
subplot(1,2,2); hold all;
title('DTV profile CVX orig under mesh refinement');
plot(idealline(:,1),idealline(:,2),'linewidth',2,'color','g');
for i=1:numel(pixels)
    plot(areas_admm{i}, perimeters_admm{i},'color', [i 0 0]/numel(pixels),'linewidth',1);
end
legend(['correct';split(num2str(pixels))]);
% exportgraphics(fh,'gammaConvTest.pdf','contenttype','vector');







