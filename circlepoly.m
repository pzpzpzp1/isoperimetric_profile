function polys = circlepoly(radii, centers, N)
    maxNumCompThreads(50);
    if nargin==0
        radii = rand(10,1);
        centers = randn(10,2);
        N=5;
    end

    theta = fliplr(linspace(0,2*pi,N+1));
    theta = theta(1:end-1);

    x = centers(:,1) + radii.*cos(theta);
    y = centers(:,2) + radii.*sin(theta);
    
    poly1 = polyshape(x(1,:),y(1,:));
    polys = repmat(poly1,numel(radii),1);
    
%     figure; hold all; axis equal;
    parfor i=1:numel(radii)
        polys(i).Vertices = [x(i,:);y(i,:)]';
    end

end