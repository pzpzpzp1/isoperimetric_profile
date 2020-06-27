function [vdb_points_pruned,vdb_edges_pruned, updown_vertices, vdb_edge_updownispoint, edgeTypes, vdb_points_pruned_multiplicity, vdb_points_pruned_radii] = process_voronoi_boundary(vdb_points, vdb_edges, vdb_edge_updownispoint, updown_vertices, points)


% g = graph(vdb_edges(:,1),vdb_edges(:,2));
% figure; axis equal; hold all;
% plot(points(:,1),points(:,2),'r.-')
% plot(g,'XData',vdb_points(:,1),'YData',vdb_points(:,2));
% vdb_points, vdb_edges, vdb_edge_updownispoint, updown_vertices, points


    % detect points outside or on boundary of polygon
    vertsInPoly = inpolygon(vdb_points(:,1),vdb_points(:,2),points(:,1),points(:,2));
    
    % detect points on boundary of polygon
    [ids, d] = knnsearch(vdb_points, points,'k',1); 
    vertsOnBoundary = zeros(size(vdb_points,1),1); vertsOnBoundary(ids)=1;
    vertsOnBoundary(ids(find(isnan(d))))=0; % false positives due to nan separators
    
    % detect concave boundary verts
    e12 = points(2:end,:)-points(1:end-1,:);
    concaveBoundaryVerts = cross([e12 e12(:,1)*0], [circshift(e12,1) e12(:,1)*0])*[0;0;1]>0;
    if sum(concaveBoundaryVerts)~=0
        [ids2, d2] = knnsearch(vdb_points, points(concaveBoundaryVerts,:)); 
        sitesOnConcaveBoundaryVerts = zeros(size(vdb_points,1),1); sitesOnConcaveBoundaryVerts(ids2)=1;
    else
        sitesOnConcaveBoundaryVerts = zeros(size(vdb_points,1),1);
    end
    % aggregate conditionals on verts
    tokeep = (vertsInPoly | vertsOnBoundary) & ~sitesOnConcaveBoundaryVerts;
    edgesToPrune = any(ismember(vdb_edges, find(~tokeep)),2);
    vdb_edges(edgesToPrune,:)=[];
    updown_vertices(edgesToPrune,:)=[];
    vdb_edge_updownispoint(edgesToPrune,:)=[];

    % trim indices.
    idx1 = 1:sum(tokeep);
    idx2 = (1:numel(tokeep))*0;
    idx2(tokeep) = idx1;
    vdb_points_pruned = vdb_points(tokeep,:);
    vdb_edges_pruned2 = sort(idx2(vdb_edges),2);
    [vdb_edges_pruned, perm] = unique(vdb_edges_pruned2,'rows');
    updown_vertices = updown_vertices(perm,:);
    vdb_edge_updownispoint = vdb_edge_updownispoint(perm,:);

    %% classify each voronoi edge.
    % project every vdb_points_pruned point onto every points edge
    [d, ~] = point_to_line_segment_distance(...
        repelem(vdb_points_pruned,size(points,1)-1,1),... 
        repmat(points(1:end-1,:), size(vdb_points_pruned,1), 1),...
        repmat(points(2:end,:), size(vdb_points_pruned,1), 1) );
    d = reshape(d, size(points,1)-1, size(vdb_points_pruned,1))';
    [vdb_points_pruned_radii] = min(d,[],2,'linear');
    vdb_points_pruned_multiplicity = sum(abs(d - vdb_points_pruned_radii)<1e-4,2);
    vdb_edges_pruned_types = sort(vdb_edge_updownispoint,2);
    vdb_edges_pruned_types_concat = vdb_edges_pruned_types(:,1)*10+vdb_edges_pruned_types(:,2);

    uniqueEdgeTypes = unique(vdb_edges_pruned_types_concat);
    
    assert(numel(uniqueEdgeTypes)<=3); % only 3 types. edge edge, vert vert, edge vert;          

    isEdgeEdge = (vdb_edges_pruned_types_concat==00);
    isVertexEdge = (vdb_edges_pruned_types_concat==01);
    isVertexVertex = (vdb_edges_pruned_types_concat==11);

    edgeTypes = [isEdgeEdge isVertexEdge isVertexVertex];

end