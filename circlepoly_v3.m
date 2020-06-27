function polys = circlepoly_v3(radii, centers, dx)
    if nargin==0
        radii = rand(10,1);
        centers = randn(10,2);
        dx = .1;
    end

    assert(numel(radii) > 0);
    polys = repmat(polyshape,numel(radii),1);
    for i=1:numel(radii)
        % more points for bigger circles.
        N = max(2*pi*radii(i)/dx, 50);
        theta = fliplr(linspace(0,2*pi,N+1));
        theta = theta(1:end-1);

        x = centers(i,1) + radii(i).*cos(theta);
        y = centers(i,2) + radii(i).*sin(theta);
    
        polys(i) = polyshape(x(1,:),y(1,:));
    end

end