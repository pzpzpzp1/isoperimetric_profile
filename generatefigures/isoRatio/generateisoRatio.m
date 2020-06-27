addpath('../../')
close all; clear all;

load pichuxys.mat;
xys(:,2) = -xys(:,2);
xys = xys - min(xys);
xys = xys./max(max(xys));
xys = xys - mean([min(xys); max(xys)]);

points1 = loadGeojson('..\..\ipvoronoi\data\All_states\Points\06_t5.geojson');
points2 = loadGeojson('..\..\ipvoronoi\data\All_states\Points\06_t2.geojson');
points1 = points1 - min(points1); points1 = points1./max(max(points1));
points1 = points1 - mean([min(points1); max(points1)]);
points2 = points2 - min(points2); points2 = points2./max(max(points2));
points2 = points2 - mean([min(points2); max(points2)]);

t = linspace(0,2*pi,1000);
r = .5+.025*sin(100*t)';
xy1 = .5*[cos(t)' sin(t)'];
xy2 = r.*[cos(t)' sin(t)'];

p1 = polyshape(xy1);
p4 = polyshape(xys);
p5 = polyshape(xy2);
p3 = polyshape(points1);
p2 = polyshape(points2);
q=@(p) 4*pi*area(p)/perimeter(p)^2;


f1=figure; axis equal; hold all; axis off; set(gcf,'color','w');
center = [0,0];
plot(center(1)+p1.Vertices([1:end,1],1),center(2)+p1.Vertices([1:end,1],2),'k'); text(center(1)-.1,center(2)+.8,num2str(q(p1)));
center = [1.2,0];
plot(center(1)+p2.Vertices([1:end,1],1),center(2)+p2.Vertices([1:end,1],2),'k'); text(center(1)-.6,center(2)+.8,num2str(q(p2)));
center = [2.4,0];
plot(center(1)+p3.Vertices([1:end,1],1),center(2)+p3.Vertices([1:end,1],2),'k'); text(center(1)-.6,center(2)+.8,num2str(q(p3)));
center = [3.6,0];
plot(center(1)+p4.Vertices([1:end,1],1),center(2)+p4.Vertices([1:end,1],2),'k'); text(center(1)-.6,center(2)+.8,num2str(q(p4)));
center = [4.8,0];
plot(center(1)+p5.Vertices([1:end,1],1),center(2)+p5.Vertices([1:end,1],2),'k'); text(center(1)-.6,center(2)+.8,num2str(q(p5)));
xlim([-.5,5.5]);
ylim([-.6,.6]);


matlab2tikz('fig1.tex','width', '.9\columnwidth','imagesAsPng',true)

