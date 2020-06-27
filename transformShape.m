function ps = transformShape(ps, scale, center, anisoScale, mockBBc)
    ps.Vertices = ps.Vertices * scale;
    ps.Vertices(:,1) = ps.Vertices(:,1)/anisoScale;
    if exist('mockBBc','var')
        % mock bounding box center
        BBc = mockBBc;
    else
        BBc = [min(ps.Vertices)+max(ps.Vertices)]/2;
    end
    ps.Vertices = ps.Vertices - BBc + center;
end