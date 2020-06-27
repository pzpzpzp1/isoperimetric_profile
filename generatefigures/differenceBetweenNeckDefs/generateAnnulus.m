clear all; close all;
t = linspace(0,pi,1000)';
xy = [cos(t) sin(t)];
fig=figure; hold all; axis equal; axis off; set(gcf,'color','w');

%% fill background
linew = 5;
dotw = 23;
neckw = 20;
x1 = [1 -.2; .95 -.25];
x2 = [.95 -.25;1 -.3];
out = [flipud([x1(:,1) x1(:,2)]);...
[[1 1]' [0 -.20]'];...
[cos(t),sin(t)];...
[[-1 -1]' [0 -.5]'];...
flipud([cos(t),-sin(t)-.5]);...
[[1 1]' [-.3 -.5]'];...
flipud([x2(:,1) x2(:,2)]);];
ps = polyshape(out);

x1 = flipud([.5 -.2; .55 -.25]);
x2 = [.55 -.25;.5 -.3];
in = [[x1(:,1),x1(:,2)];...
[.5*[1 1]' [0 -.2]'];...
[.5*cos(t),.5*sin(t)];...
[.5*[-1 -1]' [0 -.5]'];...
flipud([.5*cos(t),-.5*sin(t)-.5]);...
[.5*[1 1]' [-.3 -.5]'];...
flipud([x2(:,1),x2(:,2)])];
ps2 = polyshape(in);
psf = subtract(ps,ps2)
plot(psf,'facecolor','green','linewidth',linew)

%% fill shape
plot(cos(t),sin(t),'k','linewidth',linew)
plot(.5*cos(t),.5*sin(t),'k','linewidth',linew)
plot(.75*cos(t),.75*sin(t),'m','linewidth',linew)
plot(cos(t),-sin(t)-.5,'k','linewidth',linew)
plot(.5*cos(t),-.5*sin(t)-.5,'k','linewidth',linew)
plot(.75*cos(t),-.75*sin(t)-.5,'m','linewidth',linew)
plot([-1 -1],[0 -.5],'k','linewidth',linew)
plot(.5*[-1 -1],[0 -.5],'k','linewidth',linew)
plot(.75*[-1 -1],[0 -.5],'m','linewidth',linew)
plot([1 1],[0 -.20],'k','linewidth',linew)
plot([1 1],[-.3 -.5],'k','linewidth',linew)
plot(.5*[1 1],[0 -.2],'k','linewidth',linew)
plot(.5*[1 1],[-.3 -.5],'k','linewidth',linew)
yel = [1 0 1];
cyan = [0 1 1];
t2 = linspace(0,1,100)';
p1 = [.75 0].*t2 + (1-t2).*[.75 -.25];
cols = yel.*t2 + (1-t2).*cyan;
scatter(p1(:,1),p1(:,2),dotw,cols,'filled')
scatter(p1(:,1),-.5-p1(:,2),dotw,cols,'filled')
scatter(.75,-.25,neckw,'k','filled')
x1 = [.5 -.2; .55 -.25];
x2 = [.55 -.25;.5 -.3];
plot(x1(:,1),x1(:,2),'k','linewidth',linew)
plot(x2(:,1),x2(:,2),'k','linewidth',linew)
scatter(x1(:,1),x1(:,2),dotw,'k','filled')
scatter(x2(:,1),x2(:,2),dotw,'k','filled')
x1 = [1 -.2; .95 -.25];
x2 = [.95 -.25;1 -.3];
plot(x1(:,1),x1(:,2),'k','linewidth',linew)
plot(x2(:,1),x2(:,2),'k','linewidth',linew)
scatter(x1(:,1),x1(:,2),dotw,'k','filled')
scatter(x2(:,1),x2(:,2),dotw,'k','filled')

exportgraphics(fig,'annulus.pdf','contenttype','vector')
