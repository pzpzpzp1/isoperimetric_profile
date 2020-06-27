% compute upper bound to isoperimetric profile from opening.
% input polygon specified by points with first and last point equivalent.
% toggle for visualize.
% params = [N, dr]
% N is number of sample points on the [0,inr] interval.
% dr is the discretization size for circle polyshapes needed for interpolating
% output is perimeters and areas sorted by areas. NOT MONOTONIC.
% keypoints stores many special regions in the profile that are useful.
function [fullperimeters, fullareas, keypoints, popens, perodes] = morphOpenBound(points, Visualize, params)
    if nargin==0
        % points = loadGeojson('.\ipvoronoi\data\All_states\Points\06_t5.geojson');
        % points = userinput_points('cachedshape.lines'); % generate new shape. optionally save it.
        % points = loadCachedShape();
        points = loadPerturbedCircle(40, .025);
        points = points + [10, 40];
        Visualize = 1;
        params = [200 .04];
    end
    N = params(1); % morphological opening resolution
    dr = params(2); % morphological opening failure
    
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
    
    %% get med ax node radii ---> inr: max inscribed circle radius
    [medpoints, mededges,...
        updownpoints, updowntype,...
        edgeTypes, medpoint_degree,...
        medpoint_radii] = process_voronoi_boundary(vdb_points_raw, vdb_edges_raw, vdb_edge_updownispoint_raw, updown_vertices_raw,points);
    [inr, maxind] = max(medpoint_radii);
    inxy = medpoints(maxind,:);
    
    %% get smallest neck. quite a chunky peice of code to extract one number :(
    isEdgeEdge = edgeTypes(:,1);
    isVertexEdge = edgeTypes(:,2);
    isVertexVertex = edgeTypes(:,3);
    [ii, jj]=find(edgeTypes); [ii, perm] = sort(ii); jj = jj(perm);
    EdgesTable = table(mededges, jj, 'variablenames',{'EndNodes','EdgeTypes'});
    NodesTable = table(medpoints(:,1),medpoints(:,2),medpoints,medpoint_radii,'variablenames',{'x','y','xy','r'});
    MAG0 = graph(EdgesTable, NodesTable);
    
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
    [maxRadNeck, maxind] = max(neckr);
    [minRadNeck, mind] = min(neckr);
    minNeckXy = neckxy(mind,:);
    maxNeckXy = neckxy(maxind,:);
    numNecks = numel(neckr);
    
    if Visualize
        EdgesTable = table(unique(sort(vdb_edges_raw,2),'rows'), 'variablenames',{'EndNodes'});
        NodesTable = table(vdb_points_raw(:,1),vdb_points_raw(:,2),vdb_points_raw,'variablenames',{'x','y','xy'});
        g = graph(EdgesTable, NodesTable);
        figure; axis equal; hold all; xlim(BB(:,1)); ylim(BB(:,2)); title('full segment voronoi diagram and necks');
        plot(g,'XData',g.Nodes.x,'YData',g.Nodes.y);
        plot(points(:,1),points(:,2),'r.-');
        if numNecks~=0;
            scatter(neckxy(:,1),neckxy(:,2),'g','filled');
        end
    end
    
    %% build morphopen profile
    rads = linspace(0,inr,N);
    rads = fliplr(rads);
    popen = polyshape.empty(0);
    perode = polyshape.empty(0);
    areas = zeros(N,1);
    perimeters = zeros(N,1);
    if Visualize; figure; hold all; axis equal; plot(ps); drawnow; 
        if numNecks~=0 
            scatter(minNeckXy(1), minNeckXy(2),'k','filled'); 
        end
    end
    numErodedRegions = 1;
    newRegionEvent = [];
    for i=1:numel(rads)
        if Visualize; title([num2str(i) ' out of ' num2str(numel(rads)) ' : numErodeRegions ' num2str(numErodedRegions)]); end;
        [popen(i), perode(i)] = morphopen(ps, rads(i));
        
        areas(i) = area(popen(i));
        perimeters(i) = perimeter(popen(i));
        
        if perode(i).NumRegions > numErodedRegions
            assert(perode(i).NumRegions >= numErodedRegions + 1);
            
            % if this assert triggers it means the timestep was too big and
            % multiple regions spawned in the same timestep interval. Try
            % increasing params.nMorphOpenSamples.
            % ^ this assert message is deprecated right?
            newRegionEvent(end+1) = i;
        end
        numErodedRegions = perode(i).NumRegions;
        
        if Visualize
            try; delete(ch1); delete(ch2); catch; end;
            ch1 = plot(popen(i),'facecolor','green');
            ch2 = plot(perode(i),'facecolor','red');
            drawnow;
        end
    end
    
    %% tack on the max inscribed circle part
    startperims = (0:dr:2*pi*inr)';
    startareas = pi*(startperims/(2*pi)).^2;
    startpopen = circlepoly_v3(startperims/(2*pi),repmat(inxy,numel(startperims),1),dr);
    startperode = repmat(polyshape([0,0]),[numel(startperims),1]);
    
    %% fix the profile sampling where new region bubbles have formed.
    fixareas = [];
    fixperimeters = [];
    fixedperode = polyshape.empty(0);
    fixedpopen = polyshape.empty(0);
    for i=1:numel(newRegionEvent)
        eventInd = newRegionEvent(i);
        allregs = regions(perode(eventInd));
        
        numNewRegions = perode(eventInd).NumRegions - perode(eventInd-1).NumRegions;
        
        [~, reginds] = sort(area(allregs)); spotInds = reginds(1:numNewRegions);
        otherInds = setdiff(1:numel(allregs), spotInds);
        spotLoc = [];
        for spotInd = spotInds'
            spotLoc(end+1,:) = mean(allregs(spotInd).Vertices);
        end
        baseEroded = union(allregs(otherInds));
        baseFull = morphdilate(baseEroded, rads(eventInd));
        
        for j=0:dr:rads(eventInd)
            circ = circlepoly_v3(repmat(j,size(spotLoc,1),1), spotLoc, dr);
            E = union([baseFull; circ]);
            fixedpopen(end+1) = E;
            fixedperode(end+1) = perode(eventInd);
            fixareas(end+1) = area(E);
            fixperimeters(end+1) = perimeter(E);
        end
    end
    
    %% aggregate results
    perodes = [startperode; perode'; fixedperode'];
    popens = [startpopen; popen'; fixedpopen'];
    fullareas = [startareas; areas; fixareas'];
    fullperimeters = [startperims; perimeters; fixperimeters'];
    
    keypoints.inr = inr;
    keypoints.inxy = inxy;
    keypoints.minRadNeck = minRadNeck;
    keypoints.maxRadNeck = maxRadNeck;
    keypoints.neckxy = neckxy;
    keypoints.neckr = neckr;
    keypoints.inrArea = pi*inr^2;
    keypoints.areaMax = area(ps);
    keypoints.perimMax = perimeter(ps);
    keypoints.areaMinNeck = area(morphopen(ps, minRadNeck));
    keypoints.eventAreas = [areas(newRegionEvent-1), areas(newRegionEvent)];
    
    [fullareas, perm] = sort(fullareas);
    fullperimeters = fullperimeters(perm);
    perodes = perodes(perm);
    popens = popens(perm);
    
    %% display result
    if Visualize
        figure; hold all; xlabel('area'); ylabel('perimeter'); title('bound by opening and interoplation');
        plot(fullareas, fullperimeters,'b.-');
        xline(keypoints.inrArea,'g')
        xline(keypoints.areaMinNeck,'g')
        xline(keypoints.areaMax,'k')
        for i=1:size(keypoints.eventAreas,1)
            xline(keypoints.eventAreas(i,1),'b')
            xline(keypoints.eventAreas(i,2),'r')
        end
        
        figure; hold all; axis equal;
        plot(popens(end),'facecolor','none');
        xlim(BB(:,1)); ylim(BB(:,2))
        for i= 1:numel(perodes)
            try; delete(ch1); delete(ch2); catch; end;
            ch1 = plot(perodes(i),'facecolor','r');
            ch2 = plot(popens(i),'facecolor','g');
            drawnow;
        end
    end
end





