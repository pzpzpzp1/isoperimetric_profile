function fullenv = extractAggregateIPEnvelope(points, res, eps)
    BBcenter = [min(points)+ max(points)]/2;
    points = points - BBcenter;
    shape = polyshape(points);

    mo_areas = res.ORes.areas;
    mo_perims = res.ORes.perims;
    ma_areas = res.MARes.areas;
    ma_perims = res.MARes.perims;
    tv_areas = res.TVRes.areas;
    tv_perims = res.TVRes.perims;
    
    if numel(res.ORes.keypoints.minRadNeck)==0 || res.MARes.props.isBigNeck
        thresh = max(max(ma_perims), max(mo_perims))+1;
        % exact profile can be computed just med axis, with morph open to
        % supplement tail.
        ps1 = polyshape([[ma_areas' ma_perims']; ma_areas(end) thresh+1; 0 thresh+1;]);
        ps2 = polyshape([[mo_areas mo_perims]; mo_areas(end) thresh+1; 0 thresh+1;]);
        punion = union(ps1,ps2);
        vertsL = punion.Vertices;
        toremove = vertsL(:,2) > thresh;
        vertsL(toremove,:) = [];
        areas = vertsL(:,1);
        perims = vertsL(:,2);
        [agg_areas, reorderperm] = sort(areas);
        agg_perims = perims(reorderperm);


        vlower = [agg_areas agg_perims];
        vupper = [agg_areas agg_perims+eps];

        fullenv = polyshape([vlower; flipud(vupper)]);
        %figure; plot(agg_areas, agg_perims,'g','linewidth',2);
        return;
    end

    rn = res.ORes.keypoints.minRadNeck;
    rnArea = area(morphopen(shape, rn));

    pruneindsTV = ~(~isnan(res.TVRes.perims) & res.TVRes.perims < res.ORes.keypoints.perimMax*1.05);
    tv_areas(pruneindsTV) = [];
    tv_perims(pruneindsTV) = [];

    rt = max(res.MARes.props.maxneckrad, res.MARes.props.maxCirc2rad/2);
    inxy = mean(res.ORes.popens(2).Vertices);
    if isfield(res.ORes.keypoints, 'inxy')
        inxy = res.ORes.keypoints.inxy;
    end
    leftsidetightT = getArea1(shape, rt, inxy);
    
    thresh = max(max(ma_perims), max(mo_perims))+1;
    ps1 = polyshape([[ma_areas' ma_perims']; ma_areas(end) thresh+1; 0 thresh+1;]);
    ps2 = polyshape([[mo_areas mo_perims]; mo_areas(end) thresh+1; 0 thresh+1;]);
    punion = union(ps1,ps2);
    vertsL = punion.Vertices;
    toremove = vertsL(:,2) > thresh;
    vertsL(toremove,:) = [];
    areas = vertsL(:,1);
    perims = vertsL(:,2);
    [agg_areas, reorderperm] = sort(areas);
    agg_perims = perims(reorderperm);

    leftAreas = agg_areas(agg_areas <= leftsidetightT);
    leftPerims = agg_perims(agg_areas <= leftsidetightT);
    rightAreas = agg_areas(agg_areas >= rnArea);
    rightPerims = agg_perims(agg_areas >= rnArea);

    midAreasUp = agg_areas(agg_areas > leftsidetightT & agg_areas < rnArea);
    midPerimsUp = agg_perims(agg_areas > leftsidetightT & agg_areas < rnArea);

    midAreasDown = tv_areas(tv_areas > leftsidetightT & tv_areas < rnArea);
    midPerimsDown = tv_perims(tv_areas > leftsidetightT & tv_areas < rnArea);


    
% fix lower bound.
vertsL = [[leftAreas(end) midAreasDown rightAreas(1)];...
[leftPerims(end) midPerimsDown rightPerims(1)]]';
vertsL = [leftAreas(end) thresh+1 ;vertsL; rightAreas(1) thresh+1];
vertsU = [[leftAreas(end) midAreasUp' rightAreas(1)];...
[leftPerims(end) midPerimsUp' rightPerims(1)]]';
vertsU = [leftAreas(end) thresh+1 ;vertsU; rightAreas(1) thresh+1];

verts = union(polyshape(vertsL),polyshape(vertsU)).Vertices;
toremove = verts(:,2) > thresh;
verts(toremove,:) = [];
areas = verts(:,1);
perims = verts(:,2);
[tv_areas, reorderperm] = sort(areas);
tv_perims = perims(reorderperm);


lowerenvelope = [tv_areas, tv_perims];
upperenvelope = vertsU(2:end-1,:);


lowerenvelope(:,2) = cummax(lowerenvelope(:,2));
upperenvelope(:,2) = cummin(upperenvelope(:,2),'reverse');

%{
figure; hold all; 
plot(leftAreas, leftPerims,'b','linewidth',2);
plot(rightAreas, rightPerims,'b','linewidth',2);
plot(midAreasUp, midPerimsUp,'c','linewidth',2);
%plot(midAreasDown, midPerimsDown,'m','linewidth',2);
xline(leftsidetightT)
xline(rnArea)
% plot(tv_areas, tv_perims,'linewidth',2)
plot(lowerenvelope(:,1),lowerenvelope(:,2),'c.-','linewidth',2)
plot(upperenvelope(:,1),upperenvelope(:,2),'m.-','linewidth',2)
%}

vlower = [[leftAreas leftPerims]; lowerenvelope; [rightAreas rightPerims]];
vupper = [[leftAreas leftPerims]; upperenvelope; [rightAreas rightPerims]]; vupper(:,2) = vupper(:,2)+eps;

fullenv = polyshape([vlower; flipud(vupper)]);


end