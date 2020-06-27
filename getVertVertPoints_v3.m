function [xys, radii, isNeck, neckInds, upXYs, downXYs] = getVertVertPoints_v3(segment_point_cells, startpointc, endpointc, params, debug)
dc = params(1);
rlim = params(2);
startpoint = permute(startpointc(:,[1 2]),[2 3 1]);
endpoint = permute(endpointc(:,[1 2]),[2 3 1]);
V1 = permute(reshape(cell2mat(segment_point_cells(:,1)),2,[])',[2 3 1]);

% rotate into an easier coord system
s21 = endpoint - startpoint;
angles = atan2(s21(2,:,:),s21(1,:,:));
R = [cos(angles) sin(angles); -sin(angles) cos(angles)]; 
translationY = -multiprod(R, (startpoint+endpoint)/2); translationY = translationY(2,:,:);
translationX = -multiprod(R, V1); translationX = translationX(1,:,:);
translation = [translationX; translationY];

% rotate to easy coord system
x1 = multiprod(R, startpoint) + translation; x1 = x1(1,:,:);
x2 = multiprod(R, endpoint) + translation; x2 = x2(1,:,:);
v = multiprod(R, V1) + translation; v = v(2,:,:);
assert(all(x1(:)<=x2(:)));
isNeck = squeeze(x1<0 & x2>0);

% compute dynamic sampling
neckInds = zeros(numel(x1),1);
radii = cell(numel(x1),1);
xys = cell(numel(x1),1);
upXYs = cell(numel(x1),1);
downXYs = cell(numel(x1),1);
for i=1:numel(x1)
    x1i = x1(i);
    x2i = x2(i);
    H = abs(v(i));
    r1 = sqrt(x1i^2 + H^2);
    r2 = sqrt(x2i^2 + H^2);
    c1 = 1/r1;
    c2 = 1/r2;
    
    % choose xvals such that the difference in c is always less than or equal to dc. 
    if isNeck(i)
        cmid = 1/H;
        cvalsL = linspace(c1,cmid,1+max(ceil(abs(cmid-c1)/dc),3));
        cvalsR = linspace(cmid,c2,1+max(ceil(abs(c2-cmid)/dc),3));
        xvalsL = -sqrt(abs(1./cvalsL.^2 - H^2));
        xvalsR = sqrt(abs(1./cvalsR.^2 - H^2));
        assert(abs(xvalsL(end)) + abs(xvalsR(1)) < 1e-5);

        % based on rlim, trim the number of points.
        if numel(xvalsL)>4
            toTrimL = 1./cvalsL < rlim; toTrimL(1:2) = false; toTrimL(end) = false;
            xvalsL(toTrimL) = [];
        end
        if numel(xvalsR)>4
            toTrimR = 1./cvalsR < rlim; toTrimR(1:2) = false; toTrimR(end) = false;
            xvalsR(toTrimR) = [];
        end

        xvals = [xvalsL(1:end-1) 0 xvalsR(2:end)];
        neckInds(i) = numel(xvalsL);
    else
        cvals = linspace(c1,c2,1+max(ceil(abs(c2-c1)/dc),3));
        xvals = sqrt(abs(1./cvals.^2 - H^2)); xvals = sign(x1i)*xvals;

        % based on rlim, trim the number of points.
        if numel(xvals)>4
            toTrim = 1./cvals < rlim; toTrim(1:2) = false; toTrim(end) = false;
            xvals(toTrim) = [];
        end
    end
    
    %% store and  convert back to original coordinates
    yvals = 0*xvals;
    xyis = [xvals; yvals];
    radii{i} = sqrt(xvals.^2 + H^2);
    xys{i} = R(:,:,i)'*(xyis-translation(:,:,i));
    upXYs{i} = repmat(segment_point_cells{i,1},1,size(xys{i},2));
    downXYs{i} = repmat(segment_point_cells{i,2},1,size(xys{i},2));
    
    if debug
        p1=plot(xys{i}(1,:),xys{i}(2,:),'b.-','linewidth',2,'markersize',20);
        p2=plot(upXYs{i}(1,:),upXYs{i}(2,:),'k.-','linewidth',2,'markersize',20);
        p3=plot(downXYs{i}(1,:),downXYs{i}(2,:),'k.-','linewidth',2,'markersize',20);
        p4=plot(xys{i}(1,1),xys{i}(2,1),'r.-','linewidth',2,'markersize',20);
        p5=plot(upXYs{i}(1,1),upXYs{i}(2,1),'r.-','linewidth',2,'markersize',20);
        p6=plot(downXYs{i}(1,1),downXYs{i}(2,1),'r.-','linewidth',2,'markersize',20);
        p7=plot(xys{i}(1,end),xys{i}(2,end),'g.-','linewidth',2,'markersize',20);
        p8=plot(upXYs{i}(1,end),upXYs{i}(2,end),'g.-','linewidth',2,'markersize',20);
        p9=plot(downXYs{i}(1,end),downXYs{i}(2,end),'g.-','linewidth',2,'markersize',20);
        pause; delete(p1);    delete(p2);    delete(p3);    delete(p4);    delete(p5);    delete(p6);    delete(p7);    delete(p8);    delete(p9);
    end
end

end