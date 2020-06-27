function [fullperimeters, fullareas, xyrs, popens] = isoperim_profile_v2(points, Visualize, N, dx)
    if nargin==0
        points = loadGeojson('.\ipvoronoi\data\All_states\Points\06_t5.geojson');
        % points = userinput_points('cachedshape.lines'); % generate new shape. optionally save it.
        % points = loadCachedShape();
        Visualize = 1;
        N = 20;
        dx = .1;
    end
    
    %% recenter and bound
    BBcenter = [min(points)+ max(points)]/2;
    points = points - BBcenter;
    BB = [min(points); max(points)]*1.1;

    %% validate shape
    ps = polyshape(points,'simplify',false);
    if Visualize; figure; hold all; axis equal; plot(ps); drawnow; end;
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

    %% compute edge segment voronoi diagram (writes output to file)
    fid = fopen('temp.lines','w');
    fprintf(fid,'%f %f \n',points');
    fclose(fid);
    voronoiexe = [IPbasedir '/ipvoronoi/build/Release/ipvoronoi.exe'];
    command = strjoin({voronoiexe,'temp.lines'});
    system(command);
    delete('temp.lines')

    %% load voronoi data
    [vdb_points_raw, vdb_edges_raw, vdb_edge_updownispoint_raw, updown_vertices_raw] = loadVoronoiBoundary('ipvoronoiout.txt');

    if any(isinf(vdb_points_raw(:))) || any(isnan(vdb_points_raw(:))) ||...
            any(isinf(vdb_edges_raw(:))) || any(isnan(vdb_edges_raw(:)))
        error('Nan or inf detected!');
    end
    
    %% trim and process voronoi graph
    [medpoints, mededges,...
        updownpoints, updowntype,...
        edgeTypes, medpoint_degree,...
        medpoint_radii] = process_voronoi_boundary(vdb_points_raw, vdb_edges_raw, vdb_edge_updownispoint_raw, updown_vertices_raw,points);
    isEdgeEdge = edgeTypes(:,1);
    isVertexEdge = edgeTypes(:,2);
    isVertexVertex = edgeTypes(:,3);
    EdgeTypeColors = [.1      .1       .1;... % isEdgeEdge:     linear      %
                      .6      .6       .6;... % isVertexEdge:   quadratic   %
                      .8      .8       .8];   % isVertexVertex: hyperbolic  %

    % build coarse level medial axis graph (MAG)
    [ii, jj]=find(edgeTypes); [ii, perm] = sort(ii); jj = jj(perm);
    EdgesTable = table(mededges, jj, 'variablenames',{'EndNodes','EdgeTypes'});
    NodesTable = table(medpoints(:,1),medpoints(:,2),medpoints,medpoint_radii,'variablenames',{'x','y','xy','r'});
    MAG0 = graph(EdgesTable, NodesTable);
    
    if Visualize
        EdgesTable = table(unique(sort(vdb_edges_raw,2),'rows'), 'variablenames',{'EndNodes'});
        NodesTable = table(vdb_points_raw(:,1),vdb_points_raw(:,2),vdb_points_raw,'variablenames',{'x','y','xy'});
        g = graph(EdgesTable, NodesTable);
        figure; axis equal; hold all; xlim(BB(:,1)); ylim(BB(:,2)); title('full segment voronoi diagram');
        plot(g,'XData',g.Nodes.x,'YData',g.Nodes.y);
        plot(points(:,1),points(:,2),'r.-');
        
        figure; hold all; axis equal; xlim(BB(:,1)); ylim(BB(:,2)); title('medial axis edge types'); colorbar;
        plot(points(:,1),points(:,2),'r.-','linewidth',4);
        xA = reshape([medpoints(mededges(isEdgeEdge,1),:) medpoints(mededges(isEdgeEdge,2),:) medpoints(mededges(isEdgeEdge,1),:)*nan]',2,[])';
        plot(xA(:,1),xA(:,2),'- ','color',EdgeTypeColors(1,:),'linewidth',2);
        xA = reshape([medpoints(mededges(isVertexEdge,1),:) medpoints(mededges(isVertexEdge,2),:) medpoints(mededges(isVertexEdge,1),:)*nan]',2,[])';
        plot(xA(:,1),xA(:,2),'color',EdgeTypeColors(2,:),'linewidth',4);
        xA = reshape([medpoints(mededges(isVertexVertex,1),:) medpoints(mededges(isVertexVertex,2),:) medpoints(mededges(isVertexVertex,1),:)*nan]',2,[])';
        plot(xA(:,1),xA(:,2),'color',EdgeTypeColors(3,:),'linewidth',3);
        drawnow;
    end
    
    %% build nonlinear radius functions.
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
    
    if Visualize
        figure; hold all; axis equal; title('medial axis and radius functions'); 
        colorbar; xlim(BB(:,1)); ylim(BB(:,2)); 
        axis off; set(gcf,'color','w')
        scatter(totalXYs(1,:),totalXYs(2,:),10,totalRs(:),'filled');
        plot(points(:,1),points(:,2),'r.-');
    end
    
    totalNodeLabels = zeros(size(totalRs));
    accedges = zeros(0,2);
    accverts = zeros(0,2);
    accradii = zeros(0,1);
    for i=1:size(mededges,1)
        xye = totalXYs(:,:,i)';
        re = totalRs(:,:,i)';
        tedges = [[1:size(xye,1)-1]', [2:size(xye,1)]'];
        numCurrentNodes = size(accverts,1);
        totalNodeLabels(1,:,i) = [1:N] + numCurrentNodes;
        accedges = [accedges; tedges + numCurrentNodes];
        accverts = [accverts; xye];
        accradii = [accradii; re];
    end

    for i=1:size(medpoints,1)
        numCurrentNodes = size(accverts,1);

        newnode = medpoints(i,:);
        newnoderadii = medpoint_radii(i,:);
        accverts = [accverts; newnode];
        accradii = [accradii; newnoderadii];
        newnodelabel = numCurrentNodes + 1;

        neighborVerts = MAG0.neighbors(i);
        for j=1:numel(neighborVerts)
            mag0edge = [i neighborVerts(j)];
            [foundForward,  indF] = ismember(mag0edge, mededges,'rows');
            [foundBackward, indB] = ismember(fliplr(mag0edge), mededges,'rows');
            if foundForward
                newedge = [newnodelabel totalNodeLabels(1,1,indF)];
            else
                newedge = [newnodelabel totalNodeLabels(1,end,indB)];
            end
            accedges = [accedges; newedge];
        end
    end
    EdgesTable = table(accedges, 'variablenames',{'EndNodes'});
    NodesTable = table(accverts, accradii, 'variablenames',{'xy','r'});
    MAG1 = graph(EdgesTable, NodesTable);

    % precompute circle polys
    % candcircs = circlepoly(MAG1.Nodes.r, MAG1.Nodes.xy, circleN);
    candcircs = circlepoly_v3(MAG1.Nodes.r, MAG1.Nodes.xy, dx);
    
    %% Proceed to sample isoperimetric profile. Traverse MAG1.
    if Visualize
        figure; hold all; axis equal;
        plot(MAG1,'XData',MAG1.Nodes.xy(:,1),'YData',MAG1.Nodes.xy(:,2));
        plot(points(:,1),points(:,2),'r-','linewidth',4,'markersize',20);
    end
    [maxrad, startnode] = max(MAG1.Nodes.r);
    isIncludedNode = zeros(size(MAG1.Nodes,1),1)==1;
    isIncludedNode(startnode)=true;
    % cheegerset = circlepoly(maxrad, MAG1.Nodes.xy(startnode,:), circleN);
    cheegerset = circlepoly_v3(maxrad, MAG1.Nodes.xy(startnode,:), dx);
    MAG1A = MAG1.adjacency;
    perimeters = zeros(numel(isIncludedNode)-1,1);
    areas = zeros(numel(isIncludedNode)-1,1);
    xyrs = zeros([size(isIncludedNode,1) 3]); xyrind = 2; 
    xyrs(1,:) = [MAG1.Nodes.xy(startnode,:) maxrad];
    cheegersets = polyshape.empty(0);
    while ~all(isIncludedNode)
        try; delete(ch); delete(sh); catch; end;
        if Visualize
            cheegersets(end+1) = cheegerset;
            ch = plot(cheegerset,'facecolor','green'); 
            sh = scatter(MAG1.Nodes.xy(isIncludedNode,1),MAG1.Nodes.xy(isIncludedNode,2),20,[1 0 1],'filled');
            title([num2str(sum(isIncludedNode)) ' out of ' num2str(numel(isIncludedNode))]);
            drawnow;
        end
        currArea = area(cheegerset);
        currPerim = perimeter(cheegerset);

        % search candidates and measure their value.
        candidates = find((any(MAG1A(isIncludedNode,:),1)) & ~(isIncludedNode'));
        candidateSlope = [];
        candidateCheegers = {};
        for i=1:numel(candidates);
            candidateCircle = candcircs(candidates(i));
            candidateCheegers{i} = union(candidateCircle, cheegerset);

            candidateArea = area(candidateCheegers{i});
            candidatePerim = perimeter(candidateCheegers{i});
            candidateSlope(i) = (candidatePerim-currPerim)/(candidateArea-currArea);
        end;

        % select a candidate.
        if(any(isnan(candidateSlope)))
            selectedCandidateInd = find(isnan(candidateSlope),1);
        else
            [~, selectedCandidateInd] = min(candidateSlope);
        end

        % apply that candidate.
        selectedCandidate = candidates(selectedCandidateInd);
        cheegerset = candidateCheegers{selectedCandidateInd};
        isIncludedNode(selectedCandidate)=true;    
        perimeters(sum(isIncludedNode)-1) = perimeter(cheegerset);
        areas(sum(isIncludedNode)-1) = area(cheegerset);
        xyrs(xyrind,:) = [MAG1.Nodes.xy(candidates(selectedCandidateInd),:) MAG1.Nodes.r(candidates(selectedCandidateInd))]; xyrind = xyrind + 1;
    end
    %% augment beginning
    auga = linspace(0,areas(1),500);
    augp = 2*pi*sqrt(auga/pi);
    fullperimeters = [augp'; perimeters];
    fullareas = [auga'; areas];
    [inr, maxind] = max(medpoint_radii);
    inxy = medpoints(maxind,:);
    startpopens = circlepoly_v3(augp/(2*pi), repmat(inxy,numel(augp),1), dx);
    cheegersets = [startpopens;cheegersets'];
    popens = cheegersets;
    
    %% Visualize
    if Visualize
        figure; hold all; title('Isoperimetric Profile');
        ylabel('perimeter'); xlabel('area'); xline(areas(1),'g');
        plot(fullareas, fullperimeters,'.-'); 
    end
    
    if Visualize
        % radius function plotted out
        figure; hold all; axis equal; xlim(BB(:,1)); ylim(BB(:,2)); title('medial axis with radius-color'); colorbar;
        plot(points(:,1),points(:,2),'r.-','linewidth',4);
        %{
        xA = reshape([medpoints(mededges(isEdgeEdge,1),:) medpoints(mededges(isEdgeEdge,2),:) medpoints(mededges(isEdgeEdge,1),:)*nan]',2,[])';
        plot(xA(:,1),xA(:,2),'- ','color',EdgeTypeColors(1,:),'linewidth',2);
        xA = reshape([medpoints(mededges(isVertexEdge,1),:) medpoints(mededges(isVertexEdge,2),:) medpoints(mededges(isVertexEdge,1),:)*nan]',2,[])';
        plot(xA(:,1),xA(:,2),'color',EdgeTypeColors(2,:),'linewidth',4);
        xA = reshape([medpoints(mededges(isVertexVertex,1),:) medpoints(mededges(isVertexVertex,2),:) medpoints(mededges(isVertexVertex,1),:)*nan]',2,[])';
        plot(xA(:,1),xA(:,2),'color',EdgeTypeColors(3,:),'linewidth',3);
        %}
        scatter(xyrs(:,1), xyrs(:,2), 40, xyrs(:,3),'filled')
    end

end





