addpath('../../')
close all; clear all;
clear; close all;
load pointsopen.mat; R = [0 -1; 1 0]; points = points*R'*R'*R';

fid = fopen('temp.lines','w');
fprintf(fid,'%f %f \n',points');
fclose(fid);
voronoiexe = [IPbasedir '/ipvoronoi/build/Release/ipvoronoi.exe'];
command = strjoin({voronoiexe,'temp.lines'});
system(command);
delete('temp.lines')

    %% load voronoi data
[vdb_points_raw, vdb_edges_raw, vdb_edge_updownispoint_raw, updown_vertices_raw] = loadVoronoiBoundary('ipvoronoiout.txt');
    
%% trim and process voronoi graph
[medpoints, mededges,...
    updownpoints, updowntype,...
    edgeTypes, medpoint_degree,...
    medpoint_radii] = process_voronoi_boundary(vdb_points_raw, vdb_edges_raw, vdb_edge_updownispoint_raw, updown_vertices_raw,points);
isEdgeEdge = edgeTypes(:,1);
isVertexEdge = edgeTypes(:,2);
isVertexVertex = edgeTypes(:,3);
EdgeTypeColors = [.5      1       .5;... % isEdgeEdge:     linear      %
                  .2      .8       1;... % isVertexEdge:   quadratic   %
                  1      .7       .7];   % isVertexVertex: hyperbolic  %
              
% build coarse level medial axis graph (MAG)
[ii, jj]=find(edgeTypes); [ii, perm] = sort(ii); jj = jj(perm);
EdgesTable = table(mededges, jj, 'variablenames',{'EndNodes','EdgeTypes'});
NodesTable = table(medpoints(:,1),medpoints(:,2),medpoints,medpoint_radii,'variablenames',{'x','y','xy','r'});
MAG0 = graph(EdgesTable, NodesTable);

%% generate corase figure
f1=figure; hold all; axis equal; axis off; set(gcf,'color','white')
xA = reshape([medpoints(mededges(isEdgeEdge,1),:) medpoints(mededges(isEdgeEdge,2),:) medpoints(mededges(isEdgeEdge,1),:)*nan]',2,[])';
plot(xA(:,1),xA(:,2),'- ','color',EdgeTypeColors(1,:),'linewidth',2);
xA = reshape([medpoints(mededges(isVertexEdge,1),:) medpoints(mededges(isVertexEdge,2),:) medpoints(mededges(isVertexEdge,1),:)*nan]',2,[])';
plot(xA(:,1),xA(:,2),'color',EdgeTypeColors(2,:),'linewidth',4);
xA = reshape([medpoints(mededges(isVertexVertex,1),:) medpoints(mededges(isVertexVertex,2),:) medpoints(mededges(isVertexVertex,1),:)*nan]',2,[])';
plot(xA(:,1),xA(:,2),'color',EdgeTypeColors(3,:),'linewidth',3);
plot(points(:,1),points(:,2),'k.-','linewidth',4,'markersize',13);
scatter(medpoints(:,1),medpoints(:,2),15,[0.5977    0.2969         0],'filled')
exportgraphics(f1,'figMA_A.pdf','contenttype','vector')

%% generate fine figure
N=200;
[xys, radii, analytic] = getVertEdgePoints(updownpoints(isVertexEdge,:), medpoints(mededges(isVertexEdge,1),:), medpoints(mededges(isVertexEdge,2),:),N);
[xys2, radii2] = getVertVertPoints(updownpoints(isVertexVertex,:), medpoints(mededges(isVertexVertex,1),:), medpoints(mededges(isVertexVertex,2),:),N);
[xys3, radii3] = getEdgeEdgePoints([medpoint_radii(mededges(isEdgeEdge,1)) medpoint_radii(mededges(isEdgeEdge,2))], medpoints(mededges(isEdgeEdge,1),:), medpoints(mededges(isEdgeEdge,2),:),N);

%% combine xys, radii into one medial axis graph.
totalXYs = zeros(2,N,size(mededges,1));
totalXYs(:,:,isVertexEdge) = xys;
totalXYs(:,:,isVertexVertex) = xys2;
totalXYs(:,:,isEdgeEdge) = xys3;

totalRs = zeros(1,N,size(mededges,1));
totalRs(:,:,isVertexEdge) = radii;
totalRs(:,:,isVertexVertex) = radii2;
totalRs(:,:,isEdgeEdge) = radii3;

f2 = figure; ax1 = subplot(1,1,1);
hold all; axis equal; axis off; set(gcf,'color','w');
colormap('spring');
%colorbar; 
axis off; set(gcf,'color','w')
scatter(totalXYs(1,:),totalXYs(2,:), 20, totalRs(:),'filled');
plot(points(:,1),points(:,2),'k.-','linewidth',4,'markersize',13);
exportgraphics(f2,'figMA_B.pdf','contenttype','vector')
%cbar = fig.colorbar(ax1, [-1, 0, 1])
%matlab2tikz('MAfigB.tex','width', '.4\columnwidth')
