function [vdb_points, vdb_edges, vdb_edge_updownispoint, updown_vertices] = loadVoronoiBoundary(fname)

    if nargin == 0;
        fname = 'ipvoronoiout.txt';
    end

    fid = fopen(fname,'r');
    vdb_points = fscanf(fid,'%f');
    vdb_points = reshape(vdb_points,2,[])';
    fgetl(fid);
    vdb_edges = fscanf(fid,'%d');
    vdb_edges = reshape(vdb_edges, 4, [])';
    vdb_edge_updownispoint = vdb_edges(:,[3 4]);
    vdb_edges = vdb_edges(:,[1 2])+1;
    fgetl(fid);
    numfloatsperline = (((~vdb_edge_updownispoint)+1)*2)'; numfloatsperline=numfloatsperline(:);
    updown_vertices = fscanf(fid,'%f');
    updown_vertices = reshape(mat2cell(updown_vertices, numfloatsperline),2,[])';
    fclose(fid);

    % orient updown so that segment/point is consistent
    [vdb_edge_updownispoint, perm] = sort(vdb_edge_updownispoint,2);
    isflippedperm = sum(abs(perm - [2 1]),2)==0;
    updown_vertices(isflippedperm,:) = updown_vertices(isflippedperm,[2 1]);

end