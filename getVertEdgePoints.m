function [xys, radii, analytic, isNeck] = getVertEdgePoints(segment_point_cells, startpointc, endpointc, N)
segments = reshape(cell2mat(segment_point_cells(:,1)),4,[])';
startpoint = permute(startpointc(:,[1 2]),[2 3 1]);
endpoint = permute(endpointc(:,[1 2]),[2 3 1]);
segv2 = permute(segments(:,[1 2]),[2 3 1]);
segv1 = permute(segments(:,[3 4]),[2 3 1]);
cps = permute(reshape(cell2mat(segment_point_cells(:,2)),2,[])',[2 3 1]);

%% compute transformation needed to center locus problem onto origin.
s21 = segv2-segv1;
angles = atan2(s21(2,:,:),s21(1,:,:));
R = [cos(angles) sin(angles); -sin(angles) cos(angles)];
cv1 = cps-segv1;
projectedPoint = segv1+dot(cv1,s21,1)./vecnorm(s21,2,1).^2.*s21;
locuspoint = (projectedPoint + cps)/2;
translation = -multiprod(R, projectedPoint);
c  = multiprod(R,cps)        + translation;
flipinds = c(2,1,:) < 0; % We want to normalize the coordinate system so that c is on the positive axis. if it's negative, flip 180 and re-translate.%
pirep = ones(size(angles(flipinds)))*pi;
R(:,:,flipinds) = multiprod([cos(pirep) sin(pirep); -sin(pirep) cos(pirep)], R(:,:,flipinds));
translation = -multiprod(R, projectedPoint);

%% apply transformation (R, translation). Centers parabola at origin.
s1 = multiprod(R,segv1)      + translation;
s2 = multiprod(R,segv2)      + translation;
c  = multiprod(R,cps)        + translation;
sp = multiprod(R,startpoint) + translation;
ep = multiprod(R,endpoint)   + translation;

% parabola coefficients. y = a x^2 + k = x^2/(2c) + c/2;
cc = c(2,1,:);
a = 1./(2.*cc);
k = cc/2;

%% compute total length. 
notOriented = sp(1,1,:) >= ep(1,1,:); % 'start' is 'left' of 'end'
x1 = min([sp(1,1,:); ep(1,1,:)],[],1);
x2 = max([sp(1,1,:); ep(1,1,:)],[],1);
isNeck = x1<0 & x2>0;
L1 = x1.*sqrt(1 + x1.^2./cc.^2) + cc .* asinh(x1./cc);
L2 = x2.*sqrt(1 + x2.^2./cc.^2) + cc .* asinh(x2./cc);
len = (L2-L1)/2;

%% sample points. Eventually should use arclength sampling, but for now, just uniform N samples in x will do.
t = linspace(0,1,N);
Xs = x2.*t + x1.*(1-t);
y = a.*Xs.^2+k;
radii = y;

%% convert back to original coordinates
xy1 = [Xs;y];
xys = permute(sum(R.*permute(xy1-translation,[1 4 3 2]),1),[2 4 3 1]);
xys(:,:,notOriented) = fliplr(xys(:,:,notOriented));
radii(:,:,notOriented) = fliplr(radii(:,:,notOriented));

analytic.R = R;
analytic.translation = translation;
analytic.a = a;
analytic.k = k;
analytic.x1 = x1;
analytic.x2 = x2;
analytic.notOriented = notOriented;


%{
i=80;
figure; hold all; axis equal;
scatter(c(1,1,i),c(2,1,i),'r','filled')
scatter(s1(1,1,i),s1(2,1,i),'g','filled')
scatter(s2(1,1,i),s2(2,1,i),'g','filled')
scatter(0,0,'b','filled')
scatter(0,c(2,1,i)/2,'k','filled')
scatter(sp(1,1,i),sp(2,1,i),'k')
scatter(ep(1,1,i,1),ep(2,1,i),'k')
plot(Xs(1,:,i),y(1,:,i),'.-')
%}
end