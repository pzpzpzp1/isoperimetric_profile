% dimensions is [height width].
function [fullperimeters, fullareas, images] = getTVProfile(points, pixelsInHeight, solver, percentAreas)
    if nargin < 4
        percentAreas = linspace(0,1,20);
    end
    
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
    ps = polyshape(points);
    
    %% recenter and bound
    points = points - mean(points);
    BB = [min(points); max(points)]*1.1;
    
    %% convert to image
    totalArea = area(ps);
    totalPerimeter = perimeter(ps);
    
    yy = linspace(BB(1,2),BB(2,2),pixelsInHeight); dy = yy(2)-yy(1);
    pixelsInWidth = ceil((BB(2,1)-BB(1,1))/dy);
    xx = linspace(BB(1,1),BB(2,1),pixelsInWidth); dx = xx(2)-xx(1);
    [X,Y] = meshgrid(xx, yy);
    [binaryInd] = inpolygon(X(:), Y(:), points(:,1), points(:,2));
    image = reshape(binaryInd, pixelsInHeight, pixelsInWidth);
    % figure; axis equal; hold all; imagesc(image);
    % figure; axis equal; hold all; plot(points(:,1),points(:,2),'r','linewidth',2); 
    
    x = image;    LL = x(2:end,1:(end-1));     LR = x(2:end,2:end);     UL = x(1:(end-1),1:(end-1));     UR = x(1:(end-1),2:end); 
    diff1 = LR - LL;     diff2 = UL - LL;    diff3 = UR - UL;    diff4 = UR - LR;    v = norms([diff1(:) diff2(:) diff3(:) diff4(:)],2,2);
    L1perim = sum(v);
    
    %% compute result
    initial = image;
    images = zeros(size(image,1),size(image,2),numel(percentAreas));
    val = zeros(numel(percentAreas),1);
    for i=1:numel(percentAreas)
        percentArea = percentAreas(i);
        
        % cvx version has an improved perimeter discretization
        if strcmp(solver,'cvx') % call Mosek through CVX
            [images(:,:,i), val(i)] = isoperimetricCVX_modified(image,percentArea,[dx dy]);
        elseif strcmp(solver,'projection')
            [images(:,:,i), val(i)] = isoperimetricProjectionAlgorithm(image,percentArea,1e-5,initial);
            initial = images(:,:,i);
            val(i) = (val(i)/L1perim) * totalPerimeter;
        elseif strcmp(solver,'admm')
            [images(:,:,i), val(i)] = isoperimetricADMM(image,percentArea);
            val(i) = (val(i)/L1perim) * totalPerimeter;
        else, error('Unrecognized solver.');
        end
    end
    fullareas = percentAreas * totalArea;
    fullperimeters = val; 
end