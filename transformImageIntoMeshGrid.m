function [xyc] = transformImageIntoMeshGrid(img, scale, center, anisoScale, supposedArea)

unscale = sqrt(sum(sum(img))/supposedArea);

xs = 1:size(img,1);
ys = 1:size(img,2);
[X,Y]=meshgrid(xs,ys);

x = X(:);
y = Y(:);
col = reshape(img',[],1);
keepinds = col>.1;
x = x(keepinds);
y = y(keepinds);
xy = [y, x];
col = col(keepinds);
%figure;scatter(x,y)

xy = xy * scale / unscale;
xy(:,1) = xy(:,1)/anisoScale;
BBc = [min(xy)+max(xy)]/2;
xy = xy - BBc + center;

% col = -col;
% col = col - min(col);

xyc = [xy, col];
end