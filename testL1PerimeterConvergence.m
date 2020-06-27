clear; close all;
points = userinput_points(); 
% R = axang2rotm([0 0 1 rand*30]); points = [0 1;1 1;1 0; 0 0; 0 1]*R(1:2,1:2);

ps=polyshape(points);
exactperim = perimeter(ps);
points = points - mean(points);
BB = [min(points);max(points)] * 1.1;

figure; axis equal; hold on; plot(points(:,1),points(:,2))

pixelsInHeights = [100:100:1000 2000 3000 5000];
L1perims = zeros(numel(pixelsInHeights),1);
for i=1:numel(pixelsInHeights)
    pixelsInHeight = pixelsInHeights(i);
    yy = linspace(BB(1,2),BB(2,2),pixelsInHeight); dy = yy(2)-yy(1);
    pixelsInWidth = ceil((BB(2,1)-BB(1,1))/dy);
    xx = linspace(BB(1,1),BB(2,1),pixelsInWidth); dx = xx(2)-xx(1);
    [X,Y] = meshgrid(xx, yy);
    [binaryInd] = inpolygon(X(:), Y(:), points(:,1), points(:,2));
    image = reshape(binaryInd, pixelsInHeight, pixelsInWidth);
    
%     filter = [1 1 1; 1 1 1; 1 1 1]/9; % blur
%     image2 = conv2(image,filter);
%     image2 = image2(2:end-1,2:end-1);
    filter = ones(5); filter = filter/25; % blur
    image2 = conv2(image,filter);
    image2 = image2(3:end-2,3:end-2);
    
    x = image2;    LL = x(2:end,1:(end-1));     LR = x(2:end,2:end);     UL = x(1:(end-1),1:(end-1));     UR = x(1:(end-1),2:end); 
    diff1 = (LR - LL)/dx;     
    diff2 = (UL - LL)/dy;    
    diff3 = (UR - UL)/dx;    
    diff4 = (UR - LR)/dy;  
    
    % one version. basically just more blurring
    dfdx = (diff1 + diff3)/2;
    dfdy = (diff2 + diff4)/2;
    v1 = sqrt(dfdx(:).^2 + dfdy(:).^2);
    
    % version that was in cvx
    v2 = norms([diff1(:) diff2(:) diff3(:) diff4(:)],2,2)/sqrt(2);
    
    % no blur, 3x3 stencil
    x = image; 
    dfdx = conv2(image,[-1 -1 -1; 0 0 0; 1 1 1]'/(6*dx));
    dfdy = conv2(image,[1 1 1; 0 0 0; -1 -1 -1]/(6*dy));
    dfdx = dfdx(3:end-2);
    dfdy = dfdy(3:end-2);
    v3 = sqrt(dfdx(:).^2 + dfdy(:).^2);
    
    L1perims1(i) = sum(v1)*dx*dy;
    L1perims2(i) = sum(v2)*dx*dy;
    L1perims3(i) = sum(v3)*dx*dy;
end

figure; hold all; legend;
plot(pixelsInHeights,L1perims1); 
plot(pixelsInHeights,L1perims2); 
plot(pixelsInHeights,L1perims3); 
yline(exactperim);
ylim([0,60])
