function [inr, center] = getMaxInscribedCircle(ps, Visualize)
    if nargin<2
        Visualize = 0;
    end

    [x,y] = ps.boundary;
    points = [x,y];
    
    fid = fopen('temp.lines','w');
    fprintf(fid,'%f %f \n',points');
    fclose(fid);
    voronoiexe = [IPbasedir '/ipvoronoi/build/Release/ipvoronoi.exe'];
    command = strjoin({voronoiexe,'temp.lines'});
    system(command);
    delete('temp.lines')
    [vdb_points_raw, vdb_edges_raw, vdb_edge_updownispoint_raw, updown_vertices_raw] = loadVoronoiBoundary('ipvoronoiout.txt');

    if any(isinf(vdb_points_raw(:))) || any(isnan(vdb_points_raw(:))) ||...
            any(isinf(vdb_edges_raw(:))) || any(isnan(vdb_edges_raw(:)))
        error('Nan or inf detected!');
    end
    [medpoints, mededges,...
        updownpoints, updowntype,...
        edgeTypes, medpoint_degree,...
        medpoint_radii] = process_voronoi_boundary(vdb_points_raw, vdb_edges_raw, vdb_edge_updownispoint_raw, updown_vertices_raw,flipud(points));
    [inr,maxind]=max(medpoint_radii);
    center = medpoints(maxind,:);
    
    if Visualize
        figure; hold all; axis equal;
        plot(ps)
        patch('vertices',vdb_points_raw,'faces',vdb_edges_raw(:,[1 2 1]),'edgecolor','blue')
        patch('vertices',medpoints,'faces',mededges(:,[1 2 1]),'edgecolor','red','linewidth',2)
        plot(circlepoly_v3(inr,center,.01))
    end
        
end