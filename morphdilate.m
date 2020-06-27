function pdilate = morphdilate(ps, rad)
    if numel(rad)==0 || rad==0 || numel(ps.Vertices)==0
        pdilate = ps;
        return;
    end
    
    [x,y] = boundary(ps);
    PA = polybuffer([x,y], 'lines', rad);
    pdilate = union(ps,union(PA));

    %{
    parts = regions(ps);
    dilateparts = polyshape.empty(0,numel(parts));
    for i=1:numel(parts);
        polyj = rmholes(parts(i));
        pholes = holes(parts(i));
        for j=1:numel(pholes)
            pv = [pholes(j).Vertices; pholes(j).Vertices(1,:)];
            dilateparts(end+1) = polybuffer(pv, 'lines', rad);            
        end

        pv = [polyj.Vertices; polyj.Vertices(1,:)];
        dilateparts(end+1) = polybuffer(pv, 'lines', rad);
    end
    pdilate = union(ps,union(dilateparts));
    %}

%     figure; hold all; axis equal;
%     plot(ps);
%     plot(pdilate);
    
end