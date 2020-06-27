function fig = visualizeMedAxis(points)

%% recenter and bound
    BBcenter = [min(points)+ max(points)]/2;
    points = points - BBcenter;
    BB = [min(points); max(points)]*1.1;
    
    %% validate shape
    ps = polyshape(points,'simplify',false);
    if ~issimplified(ps)
        display('Shape errors detected. Using polyshape to simplify shape.');
        ps = polyshape(points,'simplify',true);
        points = [ps.Vertices; ps.Vertices(1,:)];
    end
    if any(isnan(points(:))) || any(isinf(points(:))) || ps.NumRegions~=1
        error('Nan or inf or split regions detected!');
    end
    totalarea = sum(cross([points(1:end-1,:) zeros(size(points,1)-1,1)],[points(2:end,:) zeros(size(points,1)-1,1)])*[0;0;1]);
    if totalarea < 0
        points = flipud(points);
    end
    
    %% get voro
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

%% get necks
[inr, maxind] = max(medpoint_radii);
inxy = medpoints(maxind,:);
isEdgeEdge = edgeTypes(:,1);
isVertexEdge = edgeTypes(:,2);
isVertexVertex = edgeTypes(:,3);
[ii, jj]=find(edgeTypes); [ii, perm] = sort(ii); jj = jj(perm);
EdgesTable = table(mededges, jj, 'variablenames',{'EndNodes','EdgeTypes'});
NodesTable = table(medpoints(:,1),medpoints(:,2),medpoints,medpoint_radii,'variablenames',{'x','y','xy','r'});
MAG0 = graph(EdgesTable, NodesTable);

%% generate corase figure
EdgesTable = table(unique(sort(vdb_edges_raw,2),'rows'), 'variablenames',{'EndNodes'});
NodesTable = table(vdb_points_raw(:,1),vdb_points_raw(:,2),vdb_points_raw,'variablenames',{'x','y','xy'});
g = graph(EdgesTable, NodesTable);
figure; axis equal; hold all; xlim(BB(:,1)); ylim(BB(:,2)); title('full segment voronoi diagram');
plot(g,'XData',g.Nodes.x,'YData',g.Nodes.y);
plot(points(:,1),points(:,2),'r.-');
        
figure; hold all; axis equal; axis off; set(gcf,'color','white')
xA = reshape([medpoints(mededges(isEdgeEdge,1),:) medpoints(mededges(isEdgeEdge,2),:) medpoints(mededges(isEdgeEdge,1),:)*nan]',2,[])';
plot(xA(:,1),xA(:,2),'- ','color',EdgeTypeColors(1,:),'linewidth',2);
xA = reshape([medpoints(mededges(isVertexEdge,1),:) medpoints(mededges(isVertexEdge,2),:) medpoints(mededges(isVertexEdge,1),:)*nan]',2,[])';
plot(xA(:,1),xA(:,2),'color',EdgeTypeColors(2,:),'linewidth',4);
xA = reshape([medpoints(mededges(isVertexVertex,1),:) medpoints(mededges(isVertexVertex,2),:) medpoints(mededges(isVertexVertex,1),:)*nan]',2,[])';
plot(xA(:,1),xA(:,2),'color',EdgeTypeColors(3,:),'linewidth',3);
plot(points(:,1),points(:,2),'k.-','linewidth',4,'markersize',13);
 
%% compute necks and fine ma
[xys, radii, analytic, isNeck, neckInds, upXYs, downXYs] = getVertEdgePoints_v3(updownpoints(isVertexEdge,:), medpoints(mededges(isVertexEdge,1),:), medpoints(mededges(isVertexEdge,2),:), [inr inr], 0);    
[xys2, radii2, isNeck2, neckInds2, upXYs2, downXYs2] = getVertVertPoints_v3(updownpoints(isVertexVertex,:), medpoints(mededges(isVertexVertex,1),:), medpoints(mededges(isVertexVertex,2),:), [inr inr], 0);
[xys3, radii3, upXYs3, downXYs3] = getEdgeEdgePoints_v3(updownpoints(isEdgeEdge,:), [medpoint_radii(mededges(isEdgeEdge,1)) medpoint_radii(mededges(isEdgeEdge,2))], medpoints(mededges(isEdgeEdge,1),:), medpoints(mededges(isEdgeEdge,2),:), [inr inr], 0);

totalXYs = cell(size(mededges,1),1);
totalXYs(isVertexEdge) = xys;
totalXYs(isVertexVertex) = xys2;
totalXYs(isEdgeEdge) = xys3;

totalNecks = zeros(size(mededges,1),1);
totalNecks(isVertexEdge) = neckInds;
totalNecks(isVertexVertex) = neckInds2;

totalRs = cell(size(mededges,1),1);
totalRs(isVertexEdge) = radii;
totalRs(isVertexVertex) = radii2;
totalRs(isEdgeEdge) = radii3;

neckxy = [];
neckr = [];
neckOi = find(isNeck); % neck outer ind
for i = 1:numel(neckOi)
    neckxy(end+1,:) = xys{neckOi(i)}(:,neckInds(neckOi(i)));
    neckr(end+1,1) = radii{neckOi(i)}(1,neckInds(neckOi(i)));
end
neckOi = find(isNeck2); % neck outer ind
for i = 1:numel(neckOi)
    neckxy(end+1,:) = xys2{neckOi(i)}(:,neckInds2(neckOi(i)));
    neckr(end+1,1) = radii2{neckOi(i)}(1,neckInds2(neckOi(i)));
end

for i=1:size(medpoints,1)
    newnoderadii = medpoint_radii(i,:);
    neighborVerts = MAG0.neighbors(i);
    neighrads = [];
    for j=1:numel(neighborVerts)
        mag0edge = [i neighborVerts(j)];
        [foundForward,  indF] = ismember(mag0edge, mededges,'rows');
        [foundBackward, indB] = ismember(fliplr(mag0edge), mededges,'rows');
        if foundForward
            neighrads(end+1) = totalRs{indF}(2);
        else
            neighrads(end+1) = totalRs{indB}(end-1);
        end
    end

    if numel(neighborVerts) > 1
        % med axis radii are unfortunately suuper close in certain non
        % generic cases. causes slight precision problems in neck
        % detection. This does mean we can't detect necks under 1e-5 so make sure the input shape is large and generic enough.
        if sum(neighrads - newnoderadii >= 1e-5) >=2 
        %if sum(neighrads >= newnoderadii) >=2 
            neckxy(end+1,:) = medpoints(i,:);
            neckr(end+1,1) = newnoderadii;
        end
    end
end
maxRadNeck = max(neckr);
[minRadNeck, mind] = min(neckr);
minNeckXy = neckxy(mind,:);
numNecks = numel(neckr);



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

fig = figure; ax1 = subplot(1,1,1);
hold all; axis equal; axis off; set(gcf,'color','w');
plot(ps,'facecolor','green')
colormap('spring');
axis off; set(gcf,'color','w')
scatter(totalXYs(1,:),totalXYs(2,:), 20, totalRs(:),'filled');
plot(points(:,1),points(:,2),'k.-','linewidth',4,'markersize',13);
if numNecks~=0
    scatter(neckxy(:,1),neckxy(:,2),40,'k','filled')
end
end