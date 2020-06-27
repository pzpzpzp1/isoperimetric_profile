function [xys, radii, analytic, isNeck, neckInds, upXYs, downXYs] = getVertEdgePoints_v3(segment_point_cells, startpointc, endpointc, params, debug)
dc = params(1);
rlim = params(2);
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
isNeck = squeeze(x1<0 & x2>0);
% L1 = x1.*sqrt(1 + x1.^2./cc.^2) + cc .* asinh(x1./cc);
% L2 = x2.*sqrt(1 + x2.^2./cc.^2) + cc .* asinh(x2./cc);
% len = (L2-L1)/2;

analytic.R = R;
analytic.translation = translation;
analytic.a = a;
analytic.k = k;
analytic.x1 = x1;
analytic.x2 = x2;
analytic.notOriented = notOriented;

neckInds = zeros(numel(x1),1);
%% sample points. 
radii = cell(numel(x1),1);
xys = cell(numel(x1),1);
upXYs = cell(numel(x1),1);
downXYs = cell(numel(x1),1);
for i=1:numel(x1)
    x1i = x1(i);
    x2i = x2(i);
    r1 = a(i)*x1i^2+k(i);
    r2 = a(i)*x2i^2+k(i);
    c1 = 1/r1;
    c2 = 1/r2;
    
    % choose xvals such that the difference in c is always less than or equal to dc. 
    if isNeck(i)
        cmid = 1/k(i);
        cvalsL = linspace(c1,cmid,1+max(ceil(abs(cmid-c1)/dc),3));
        cvalsR = linspace(cmid,c2,1+max(ceil(abs(c2-cmid)/dc),3));
        %assert(all((1./cvalsL - k(i)) >= 0)); % some floating precision makes this sometimes false.
        xvalsL = sqrt(abs(1./cvalsL - k(i))/a(i)); xvalsL = -xvalsL;
        %assert(all((1./cvalsR - k(i)) >= 0)); % some floating precision makes this sometimes false.
        xvalsR = sqrt(abs(1./cvalsR - k(i))/a(i)); 
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
        
        xvals = [xvalsL xvalsR(2:end)];
        neckInds(i) = numel(xvalsL);
    else
        cvals = linspace(c1,c2,1+max(ceil(abs(c2-c1)/dc),3));
        assert(all((1./cvals - k(i)) >= 0)); % should always be non-negative
        xvals = sqrt((1./cvals-k(i))/a(i)); xvals = sign(x1i)*xvals;

        % based on rlim, trim the number of points.
        if numel(xvals)>4
            toTrim = 1./cvals < rlim; toTrim(1:2) = false; toTrim(end) = false;
            xvals(toTrim) = [];
        end
    end
    
    %% store and  convert back to original coordinates
    yvals = a(i)*xvals.^2+k(i);
    xyis = [xvals; yvals];
    upXYis = [xvals; 0*yvals];
    downXYis = c(:,1,i);
    radii{i} = yvals;
    xys{i} = R(:,:,i)'*(xyis-translation(:,:,i));
    upXYs{i} = R(:,:,i)'*(upXYis-translation(:,:,i));
    downXYs{i} = repmat(R(:,:,i)'*(downXYis-translation(:,:,i)),1,size(upXYs{i},2));
    
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
    
    if notOriented(i)
        radii{i} = fliplr(radii{i});
        xys{i} = fliplr(xys{i});
        upXYs{i} = fliplr(upXYs{i});
        downXYs{i} = fliplr(downXYs{i});
    end
end

end