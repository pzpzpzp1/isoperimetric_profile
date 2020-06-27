clear all; close all;

%% build shapes
points = loadPerturbedCircle(40, .025);
pertCircPs = polyshape(points);
newrad = sqrt(area(pertCircPs)/pi);
CircPs = circlepoly_v3(newrad,[0,0],1000);

%% compute profiles
ppv = [pertCircPs.Vertices;pertCircPs.Vertices(1,:)];
[fullperimeters, fullareas, keypoints] = morphOpenBound(ppv, 1, [100, .04]);

circperims = linspace(0,2*pi*newrad,1000);
circareas = (circperims/(2*pi)).^2*pi;

%% get convex lower envelope

areasL = linspace(0, keypoints.inrArea,10);
perimsL = linspace(0, 2*pi*keypoints.inr,10);
[~,ind] = min(abs(fullareas - keypoints.inrArea));
areasR = fullareas(ind:end,:);
perimsR = fullperimeters(ind:end,:);
cvxLAreas = [areasL areasR'];
cvxLPerims = [perimsL perimsR'];


%% visualize results
fh=figure; hold all; set(gcf,'color','w');
xlabel('Area'); ylabel('Perimeter'); title('Isoperimetric Profile of Circle and Perturbed Circle');
plot([-10 0], [-10 0], 'r','linewidth',2);
plot([-10 0], [-10 0], 'g','linewidth',2);
plot([-10 0], [-10 0], 'b','linewidth',2);
xlim([0, area(pertCircPs)]);
ylim([0, perimeter(pertCircPs)]);
fh.OuterPosition=([1        1          800        350.4]);
p1 = plot(cvxLAreas, cvxLPerims-.02, 'b','linewidth',2);
p2 = plot(fullareas, fullperimeters, 'r','linewidth',2);
p3 = plot(circareas, circperims-.02, 'g','linewidth',2);

% Because the circles are on different scales
% from perim versus area, had to rescale them so the plot shows them as
% circles and not ellipses. The scaling only works if the figure is docked
% on the bottom left. some scaling parameter i must have missed.
% obviously , the circles are not on the same scale as
% the profile.
rscale = 1.5;
centerN = [.6, 1.5];
centerP = [.2, 3.3];
axvals = axis;
axisScale = (axvals(4)-axvals(3))/(axvals(2)-axvals(1));
axisScale = axisScale/((fh.OuterPosition(4))/(fh.OuterPosition(3)));
pertCircPsV = polyshape(points);
CircPsV = circlepoly_v3(newrad,[0,0],1000);
pertCircPsV.Vertices = rscale*pertCircPsV.Vertices;
CircPsV.Vertices = rscale*CircPsV.Vertices;
pertCircPsV.Vertices(:,1) = pertCircPsV.Vertices(:,1)/axisScale;
CircPsV.Vertices(:,1) = CircPsV.Vertices(:,1)/axisScale;
pertCircPsV.Vertices = pertCircPsV.Vertices + centerP;
CircPsV.Vertices = CircPsV.Vertices + centerN;
plot(CircPsV,'facecolor','green')
plot(pertCircPsV,'facecolor','red')
fh.OuterPosition=([1        1          600        350.4]);
lgd = legend({'perturbed circle','circle','cvx lower p.c.'});
lgd.Position = [0.6215    0.6767    0.2223    0.1940];
exportgraphics(fh,'pertVFull.pdf','contenttype','vector');



