function perode = morpherode(ps, rad)
    if numel(rad)==0 || rad==0 || numel(ps.Vertices)==0
        perode = ps;
        return;
    end
    
    [x,y] = boundary(ps);
    PA = polybuffer([x,y], 'lines', rad);
    perode = subtract(ps,union(PA));

    %{
    parts = regions(ps);
    erodeparts = polyshape.empty(0,numel(parts));
    for i=1:numel(parts);
        polyj = rmholes(parts(i));

        pholes = holes(parts(i));
        for j=1:numel(pholes)
            ph = [pholes(j).Vertices; pholes(j).Vertices(1,:)];
            erodeparts(end+1) = polybuffer(ph, 'lines', rad);
        end
        
        pv = [polyj.Vertices; polyj.Vertices(1,:)];
        erodeparts(end+1) = polybuffer(pv, 'lines', rad);
    end
    perode = subtract(ps,union(erodeparts));
    %}

%     figure; hold all; axis equal;
%     plot(ps);
%     plot(perode);
    
end