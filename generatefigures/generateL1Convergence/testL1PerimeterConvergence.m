clear; close all;
load pointsopen.mat;
% points = userinput_points();
R = axang2rotm([0 0 1 pi/2]); 
points = points*R(1:2,1:2);
ps = polyshape(points);
exactperim = perimeter(ps);
points = points - mean(points);
BB = [min(points);max(points)] * 1.1;

pixelsInHeights = [10:10:90 100:100:1000 2000];
L1perims = zeros(numel(pixelsInHeights),1);
for i=1:numel(pixelsInHeights)
    pixelsInHeight = pixelsInHeights(i);
    yy = linspace(BB(1,2),BB(2,2),pixelsInHeight); dy = yy(2)-yy(1);
    pixelsInWidth = ceil((BB(2,1)-BB(1,1))/dy);
    xx = linspace(BB(1,1),BB(2,1),pixelsInWidth); dx = xx(2)-xx(1);
    [X,Y] = meshgrid(xx, yy);
    [binaryInd] = inpolygon(X(:), Y(:), points(:,1), points(:,2));
    image = reshape(binaryInd, pixelsInHeight, pixelsInWidth);
    
    % All methods perform reasonably when the image is blurred out.
%     filter = ones(5); filter = filter/25; % blur
%     image = conv2(image,filter);
%     image = image(3:end-2,3:end-2);
    
    
    if i==1 
        %f1 = figure; hold all; axis equal; imagesc(image); imshow(~image);
        %xlim([0,size(image,2)]); ylim([0,size(image,1)]); axis off; set(gcf,'color','w');
        %saveas(image,'lowRes.png', 'png')
        img = imresize(~image,pixelsInHeights(end)/pixelsInHeights(1),'method','box');
        inds = find(~all(img==1));
        img2 = img(:,min(inds):max(inds));
        inds2 = find(~all(img'==1));
        img3 = img2(min(inds2):max(inds2),:);
        imwrite(img3,'lowRes.png');
    elseif i==numel(pixelsInHeights)
        %f2 = figure; hold all; axis equal; imagesc(image); imshow(~image); 
        %xlim([0,size(image,2)]); ylim([0,size(image,1)]); axis off; set(gcf,'color','w');
        img = ~image;
        inds = find(~all(img==1));
        img2 = img(:,min(inds):max(inds));
        inds2 = find(~all(img'==1));
        img3 = img2(min(inds2):max(inds2),:);
        imwrite(img3,'highRes.png');
    end
    
    % version that was in cvx
    x = image;    LL = x(2:end,1:(end-1));     LR = x(2:end,2:end);     UL = x(1:(end-1),1:(end-1));     UR = x(1:(end-1),2:end); 
    diff1 = (LR - LL)/dx;     
    diff2 = (UL - LL)/dy;    
    diff3 = (UR - UL)/dx;    
    diff4 = (UR - LR)/dy;  
    v2 = norms([diff1(:) diff2(:) diff3(:) diff4(:)],2,2)/sqrt(2); % triple checked the sqrt(2) is the right scaling
    
    % 3x3 stencil
    dfdx = conv2(image,[-1 -1 -1; 0 0 0; 1 1 1]'/(6*dx));
    dfdy = conv2(image,[1 1 1; 0 0 0; -1 -1 -1]/(6*dy));
    dfdx = dfdx(3:end-2);
    dfdy = dfdy(3:end-2);
    v3 = sqrt(dfdx(:).^2 + dfdy(:).^2);
    
    % 5x5 stencil
    dfdxFilter = repmat([1/12 	-2/3 	0 	2/3 	-1/12]/(5*dx),5,1);
    dfdyFilter = flipud(repmat([1/12 	-2/3 	0 	2/3 	-1/12]/(5*dy),5,1)');
    dfdx = conv2(image,dfdxFilter);
    dfdy = conv2(image,dfdyFilter);
    v1 = sqrt(dfdx(:).^2 + dfdy(:).^2);
    
    L1perims1(i) = sum(v1)*dx*dy;
    L1perims2(i) = sum(v2)*dx*dy;
    L1perims3(i) = sum(v3)*dx*dy;
    L1perims4(i) = sum(sqrt(diff1(:).^2+diff2(:).^2))*dx*dy;
end

fh=figure; hold all; 
fh.OuterPosition = [1            1        654.4        360.8];
yline(exactperim,'g','linewidth',2);
plot(pixelsInHeights, L1perims2,'linewidth',2);  % old cvx version
plot(pixelsInHeights, L1perims3,'linewidth',2); % new cvx version 3x3
plot(pixelsInHeights, L1perims1,'linewidth',2);  % old cvx version 5x5
legend({'Exact perimeter','TV-1','3x3 Gradient','5x5 Gradient'});
title('Convergence of discrete TV to perimeter');
xlabel('Horizontal Pixel Resolution');
ylabel('Discrete TV');


fh=figure; hold all; 
fh.OuterPosition = [1            1        654.4        360.8];
yline(exactperim,'g','linewidth',2);
plot(pixelsInHeights, L1perims4,'linewidth',2);  % chambolle survey 'correct'
plot(pixelsInHeights, L1perims2,'linewidth',2);  % old cvx version
plot(pixelsInHeights, L1perims3,'linewidth',2); % new cvx version 3x3
plot(pixelsInHeights, L1perims1,'linewidth',2);  % old cvx version 5x5
legend({'Exact perimeter','TV-0','TV-1','3x3 Gradient','5x5 Gradient'});
title('Convergence of discrete TV to perimeter');
xlabel('Horizontal Pixel Resolution');
ylabel('Discrete TV');





