function [xys, radii, upXYs, downXYs] = getEdgeEdgePoints_v3(updowns, radii_start_and_end, startpointc, endpointc, params, debug)
startpoint = permute(startpointc(:,[1 2]),[2 3 1]);
endpoint = permute(endpointc(:,[1 2]),[2 3 1]);
dc = params(1);
rlim = params(2);

%% sample points. 
radii = cell(size(radii_start_and_end,1),1);
xys = cell(size(radii_start_and_end,1),1);
upXYs = cell(size(radii_start_and_end,1),1);
downXYs = cell(size(radii_start_and_end,1),1);
for i=1:size(radii_start_and_end,1)
    xy1 = startpoint(:,1,i);
    xy2 = endpoint(:,1,i);
    r1 = radii_start_and_end(i,1);
    r2 = radii_start_and_end(i,2);
    
    if r1<rlim && r2<rlim
        % features are too small. numerically too expensive to capture
        % them. turn up the resolution here to forcefully capture it. it will cost a lot of time.
        display('Entire MA edge has features too small to capture!');
        tvals = [0 .25 .5 .75 1];
    elseif r1>rlim && r2>rlim
        % all good. all features are above the limit.
        c1 = 1/r1;
        c2 = 1/r2;
        tvals = linspace(0, 1, 1+max(ceil(abs(c1-c2)/dc),3));
    else
        tlim = (rlim-r1)/(r2-r1);
        if r1 < r2
            cvals = [linspace(1/rlim, 1/r2, 1+max(ceil(abs(1/rlim - 1/r2)/dc),3))];
            rvals = [r1 1./cvals];
            % rvals = [1./cvals]; % rlim so exclude small feature
        else
            cvals = [linspace(1/r1, 1/rlim, 1+max(ceil(abs(1/r1 - 1/rlim)/dc),3))];
            rvals = [1./cvals r2];
            % rvals = [1./cvals]; % rlim so exclude small feature
        end
        
        tvals = (rvals-r1)/(r2-r1);      
    end
    radii{i} = r1*(1-tvals) + r2*tvals;
    xys{i} = xy1*(1-tvals) + xy2*tvals;
    
    % project xys to up-down lines. 
    s1 = updowns{i,1}(1:2); s2 = updowns{i,1}(3:4);
    e12 = s2-s1; e12 = e12/norm(e12); x12 = xys{i} - s1; proj1 = (e12'*x12).*e12 + s1;
    s1 = updowns{i,2}(1:2); s2 = updowns{i,2}(3:4);
    e12 = s2-s1; e12 = e12/norm(e12); x12 = xys{i} - s1; proj2 = (e12'*x12).*e12 + s1;
    upXYs{i} = proj1;
    downXYs{i} = proj2;
    
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