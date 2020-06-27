function [fullperimeters, fullareas, cheegersets, props] = isoperim_profile_v3(points, Visualize, params, debug)

    if nargin==0
        points = loadGeojson('.\ipvoronoi\data\All_states\Points\06_t5.geojson');
        % points = userinput_points('cachedshape.lines'); % generate new shape. optionally save it.
        % points = loadCachedShape();
        Visualize = 1;
        params.dc = .1;
        params.rlim = .01;
        params.dx = .1;
        params.postprocess = 0;
        params.postprocessclosingradius = nan;
    end
    dc = params.dc; % change in candidacy
    rlim = params.rlim; % minimum radius
    dx = params.dx; % edge size of circle approx
    
    %% recenter and bound
    BBcenter = (min(points)+ max(points))/2;
    points = points - BBcenter;
    BB = [min(points); max(points)]*1.1;

    %% validate shape
    ps = polyshape(points,'simplify',false);
    if Visualize; figure; hold all; axis equal; plot(ps); drawnow; end;
    if ~issimplified(ps)
        disp('Shape errors detected. Using polyshape to simplify shape.');
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
    if debug; f1=figure; hold all; axis equal; plot(points(:,1),points(:,2),'r-','linewidth',1); end;
    [xys, radii, analytic, isNeck, neckInds, upXYs, downXYs] = getVertEdgePoints_v3(updownpoints(isVertexEdge,:), medpoints(mededges(isVertexEdge,1),:), medpoints(mededges(isVertexEdge,2),:), [dc rlim], debug);    
    [xys2, radii2, isNeck2, neckInds2, upXYs2, downXYs2] = getVertVertPoints_v3(updownpoints(isVertexVertex,:), medpoints(mededges(isVertexVertex,1),:), medpoints(mededges(isVertexVertex,2),:), [dc rlim], debug);
    [xys3, radii3, upXYs3, downXYs3] = getEdgeEdgePoints_v3(updownpoints(isEdgeEdge,:), [medpoint_radii(mededges(isEdgeEdge,1)) medpoint_radii(mededges(isEdgeEdge,2))], medpoints(mededges(isEdgeEdge,1),:), medpoints(mededges(isEdgeEdge,2),:), [dc rlim], debug);
    if debug; close(f1); end;
    
    %% combine xys, radii into one medial axis graph.
    totalUPs = cell(size(mededges,1),1);
    totalUPs(isVertexEdge) = upXYs;
    totalUPs(isVertexVertex) = upXYs2;
    totalUPs(isEdgeEdge) = upXYs3;
    
    totalDOWNs = cell(size(mededges,1),1);
    totalDOWNs(isVertexEdge) = downXYs;
    totalDOWNs(isVertexVertex) = downXYs2;
    totalDOWNs(isEdgeEdge) = downXYs3;
    
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
    
    if Visualize && false
        figure; hold all; axis equal; title('medial axis and radius functions'); 
        colorbar; xlim(BB(:,1)); ylim(BB(:,2)); 
        axis off; set(gcf,'color','w')
        for i=1:numel(totalXYs)
            scatter(totalXYs{i}(1,:),totalXYs{i}(2,:),10,totalRs{i}(:),'filled');
            if totalNecks(i)~=0
                scatter(totalXYs{i}(1,totalNecks(i)),totalXYs{i}(2,totalNecks(i)),30,'k','filled');
            end
        end
        plot(points(:,1),points(:,2),'r.-');
    end
    
    totalNodeLabels = cell(size(totalRs,1),1);
    accedges = zeros(0,2); % accumulated vals
    accverts = zeros(0,2);
    accradii = zeros(0,1);
    neckVerts = [];
    vert2mededge = zeros(0,2);
    edge2mededge = zeros(0,1);
    edge2poly = polyshape.empty(0,1);
    for i=1:size(mededges,1)
        % prune endpoints because they will be added later
        xye = totalXYs{i}(:,:)'; xye = xye(2:end-1,:);
        re = totalRs{i}(:,:)'; re = re(2:end-1,:);
        tedges = [[1:size(xye,1)-1]', [2:size(xye,1)]'];
        
        numCurrentNodes = size(accverts,1);
        totalNodeLabels{i}(1,:) = [1:size(xye,1)] + numCurrentNodes;
        accedges = [accedges; tedges + numCurrentNodes];
        edge2mededge = [edge2mededge; repmat(i,size(tedges,1),1)];
        accverts = [accverts; xye];
        accradii = [accradii; re];
        
%         figure; hold all; axis equal; title('test edge polygons'); 
%         colorbar; xlim(BB(:,1)); ylim(BB(:,2)); 
%         axis off; set(gcf,'color','w')
%         plot(points(:,1),points(:,2),'r-','linewidth',1,'markersize',20);
        for j = 1:size(tedges,1)
            edge2poly(end+1,1) = polyshape([totalUPs{i}(:,tedges(j,:)+1), fliplr(totalDOWNs{i}(:,tedges(j,:)+1))]');
            
%             p2 = plot(edge2poly(end),'facecolor','green');
%             vs = accverts(accedges(numel(edge2poly),:),:);
%             p1 = plot(vs(:,1),vs(:,2),'k.-','linewidth',1);
%             pause;
%             delete(p1); delete(p2);
        end
        
        vert2mededge(numCurrentNodes+(1:size(xye,1)),1) = i; 
        vert2mededge(numCurrentNodes+(1:size(xye,1)),2) = (1:size(xye,1))+1;
        if totalNecks(i)
            % keep track of necks. -1 is for the pruned endpoints
            neckVerts(end+1) = numCurrentNodes + totalNecks(i) - 1;
        end
    end

    isMedpointVert = false(size(medpoints,1) + size(accverts,1), 1);
    for i=1:size(medpoints,1)
        numCurrentNodes = size(accverts,1);

        newnode = medpoints(i,:);
        newnoderadii = medpoint_radii(i,:);
        accverts = [accverts; newnode];
        accradii = [accradii; newnoderadii];
        newnodelabel = numCurrentNodes + 1;
        isMedpointVert(newnodelabel) = true;

        neighborVerts = MAG0.neighbors(i);
        neighrads = [];
        
%         figure; hold all; axis equal; title('test edge polygons'); 
%         colorbar; xlim(BB(:,1)); ylim(BB(:,2)); 
%         axis off; set(gcf,'color','w')
%         plot(points(:,1),points(:,2),'r-','linewidth',1,'markersize',20);
        for j=1:numel(neighborVerts)
            mag0edge = [i neighborVerts(j)];
            [foundForward,  indF] = ismember(mag0edge, mededges,'rows');
            [foundBackward, indB] = ismember(fliplr(mag0edge), mededges,'rows');
            if foundForward
                newedge = [newnodelabel totalNodeLabels{indF}(1,1)];
                neighrads(end+1) = accradii(totalNodeLabels{indF}(1,1));
                polyverts = [totalUPs{indF}(:,1:2) fliplr(totalDOWNs{indF}(:,1:2))]';
            else
                newedge = [newnodelabel totalNodeLabels{indB}(1,end)];
                neighrads(end+1) = accradii(totalNodeLabels{indB}(1,end));
                polyverts = [totalUPs{indB}(:,end-1:end) fliplr(totalDOWNs{indB}(:,end-1:end))]';
            end
            accedges = [accedges; newedge];
            edge2mededge = [edge2mededge; indB+indF]; % one of them is always 0 so we can add them.
            edge2poly(end+1) = polyshape(polyverts);
            
%             p2 = plot(edge2poly(end),'facecolor','green');
%             vs = accverts(accedges(numel(edge2poly),:),:);
%             p1 = plot(vs(:,1),vs(:,2),'g.-','linewidth',1);
%             pause;
%             delete(p1); delete(p2);
        end
        
        if numel(neighborVerts) > 1
            % detect neck
            if sum(neighrads > newnoderadii) >=2 
                neckVerts(end+1) = newnodelabel;
                display('Found neck at coarse MA-node');
            end
        end
    end
    EdgesTable = table(accedges, edge2poly,edge2mededge, 'variablenames',{'EndNodes','Poly','MAedge'});
    NodesTable = table(accverts, accradii, 'variablenames',{'xy','r'});
    MAG1 = graph(EdgesTable, NodesTable);

    if Visualize
        figure; hold all; axis equal; title('MAG1: medial axis and radius functions'); 
        colorbar; xlim(BB(:,1)); ylim(BB(:,2)); 
        axis off; set(gcf,'color','w')
        plot(MAG1,'XData',MAG1.Nodes.xy(:,1),'YData',MAG1.Nodes.xy(:,2));
        plot(points(:,1),points(:,2),'r-','linewidth',4,'markersize',20);
        scatter(accverts(:,1),accverts(:,2),10,accradii,'filled');
        scatter(accverts(neckVerts,1),accverts(neckVerts,2),30,'k','filled');
    end
    
    if debug
        figure; hold all; axis equal; title('test edge polygons'); 
        colorbar; xlim(BB(:,1)); ylim(BB(:,2)); 
        axis off; set(gcf,'color','w')
        plot(points(:,1),points(:,2),'r-','linewidth',1,'markersize',20);
        
        plot(points(:,1),points(:,2),'r.-','linewidth',4);
        xA = reshape([medpoints(mededges(isEdgeEdge,1),:) medpoints(mededges(isEdgeEdge,2),:) medpoints(mededges(isEdgeEdge,1),:)*nan]',2,[])';
        plot(xA(:,1),xA(:,2),'- ','color',EdgeTypeColors(1,:),'linewidth',2);
        xA = reshape([medpoints(mededges(isVertexEdge,1),:) medpoints(mededges(isVertexEdge,2),:) medpoints(mededges(isVertexEdge,1),:)*nan]',2,[])';
        plot(xA(:,1),xA(:,2),'color',EdgeTypeColors(2,:),'linewidth',4);
        xA = reshape([medpoints(mededges(isVertexVertex,1),:) medpoints(mededges(isVertexVertex,2),:) medpoints(mededges(isVertexVertex,1),:)*nan]',2,[])';
        plot(xA(:,1),xA(:,2),'color',EdgeTypeColors(3,:),'linewidth',3);
        drawnow;
        
        for i=1:size(MAG1.Edges.EndNodes,1)
            type = find(edgeTypes(MAG1.Edges.MAedge(i),:));
            title([num2str(i) ' type:' num2str(type)])
            vs = MAG1.Nodes.xy(MAG1.Edges.EndNodes(i,:),:);
            if type == 3
                % should have no polyshape. vert-vert case.
                p1 = plot(vs(:,1),vs(:,2),'r.-','linewidth',1);
            else
                p1 = plot(vs(:,1),vs(:,2),'g.-','linewidth',1);
            end
            p2 = plot(MAG1.Edges.Poly(i),'facecolor','green');
            
            pause;
            delete(p1); delete(p2);
        end
    end
    
    %% compute if we are a 'big neck' shape
    minneckrad = min(accradii(neckVerts));
    maxneckrad = max(accradii(neckVerts));
    [maxrad, startnode] = max(MAG1.Nodes.r);
    maxCirc1 = circlepoly_v3(maxrad, MAG1.Nodes.xy(startnode,:), dx);
    remainder = subtract(ps,maxCirc1);
    remainderRegions = regions(remainder);
    maxCirc2rad = -1;
    maxCirc2Center = [nan nan];
    for j=1:numel(remainderRegions)
        [radius, center] = getMaxInscribedCircle(remainderRegions(j),Visualize);
        if radius > maxCirc2rad
            maxCirc2rad = radius;
            maxCirc2Center = center;
        end
    end
    props.minneckrad = minneckrad;
    props.maxneckrad = maxneckrad;
    props.maxCirc2rad = maxCirc2rad;
    props.maxCirc2Center = maxCirc2Center;
    props.isBigNeck = props.minneckrad > props.maxCirc2rad/2;
    
    %% precompute circle polys
    medpointCircles = circlepoly_v3(medpoint_radii, medpoints, dx);
    all_circs = circlepoly_v3(MAG1.Nodes.r, MAG1.Nodes.xy, dx);
                
    %% Proceed to sample isoperimetric profile. Traverse MAG1.
    if Visualize
        figure; hold all; axis equal;
        plot(MAG1,'XData',MAG1.Nodes.xy(:,1),'YData',MAG1.Nodes.xy(:,2));
        plot(points(:,1),points(:,2),'r-','linewidth',4,'markersize',20);
    end
    [maxrad, startnode] = max(MAG1.Nodes.r);
    isIncludedNode = false(size(MAG1.Nodes,1),1);
    isIncludedNode(startnode)=true;
    isIncludedEdge = false(size(MAG1.Edges,1),1);
    cheegersets(1) = circlepoly_v3(maxrad, MAG1.Nodes.xy(startnode,:), dx);
    MAG1A = MAG1.adjacency;
    perimeters = perimeter(cheegersets(1));
    areas = area(cheegersets(1));
    while ~all(isIncludedEdge) 
        currArea = area(cheegersets(end));
        currPerim = perimeter(cheegersets(end));

        % search through candidates and measure their value.
        candidateSlope = [];
        candidateAddedEdge = [];
        candidateAddedNode = [];
        numCandEdges = [];
        candI = [];
        for i=find(isIncludedNode)'
            [adjEdges, adjNodes] = outedges(MAG1,i);
            candidateEdges = find(~ismember(adjEdges, find(isIncludedEdge)));
            for j=candidateEdges(:)'
                numCandEdges(end+1) = numel(candidateEdges);
                candI(end+1) = i;
                candidateAddedEdge(end+1) = adjEdges(j);
                candidateAddedNode(end+1) = adjNodes(j);
                candSlope = 1/MAG1.Nodes.r(candidateAddedNode(end));
                candidateSlope(end+1) = candSlope;
            end
        end
        % select a candidate.
        [~, selectedCandidateInd] = min(candidateSlope);
        

        % build candidate
        edgeInd = candidateAddedEdge(selectedCandidateInd);
        nodeInd = candidateAddedNode(selectedCandidateInd);
        isIncludedNode(nodeInd) = 1;
        if numCandEdges(selectedCandidateInd)==1
            % no more outgoing edges so node is not needed anymore.
            isIncludedNode(candI(selectedCandidateInd)) = 0;
        end
        isIncludedEdge(edgeInd) = 1;
        % consider dynamically updating instead of recomputing
        coveredMedpointVerts = intersect(unique(MAG1.Edges.EndNodes(isIncludedEdge,:)),find(isMedpointVert));
        vertCircles = unique([find(isIncludedNode); coveredMedpointVerts]);
        % shapes2union needs to have empty polyshapes removed first because matlab union can have buggy behavior otherwise :(.
        shapes2union = [MAG1.Edges.Poly(isIncludedEdge); all_circs(vertCircles)];
        selectedCandidateShape = union(shapes2union(area(shapes2union)~=0));
        
        % postprocess shape
        if isfield(params,'postprocess') && params.postprocess
            assert(isfield(params,'postprocessclosingradius'));

            selectedCandidateShape = postProcessCandidatePoly(ps, selectedCandidateShape, params.postprocessclosingradius);            
        end
        
        % store results of that candidate
        cheegersets(end+1) = selectedCandidateShape;
        perimeters(end+1) = perimeter(cheegersets(end));
        areas(end+1) = area(cheegersets(end));
        
        % visualize end state of propagation
        try; delete(ch); delete(sh); catch; end;
        if Visualize
            ch = plot(cheegersets(end),'facecolor','green'); 
            sh = scatter(MAG1.Nodes.xy(isIncludedNode,1),MAG1.Nodes.xy(isIncludedNode,2),20,[1 0 1],'filled');
            title([num2str(sum(isIncludedEdge)) ' out of ' num2str(numel(isIncludedEdge))]);
            drawnow;
        end
    end
    
    
    % augment beginning. uniform sampling on perimeter produces better approximation to sqrt than uniform sampling on area.
    pmax = 2*pi*sqrt(areas(1)/pi); 
    augp = 0:.2:pmax;
    auga = (augp./(2*pi)).^2*pi; 
    fullperimeters = [augp, perimeters];
    fullareas = [auga, areas];
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
        figure; hold all; axis equal; xlim(BB(:,1)); ylim(BB(:,2)); title('medial axis with radius-color'); 
        plot(points(:,1),points(:,2),'r.-','linewidth',4);
        for i=1:numel(cheegersets)
            try; delete(p1);catch;end;
            p1 = plot(cheegersets(i),'facecolor','green');
            drawnow;
        end
    end

end





