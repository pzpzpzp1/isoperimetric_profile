function fh=visualizeProfile(points, res, params, visparams)

    if nargin < 4
        visparams.nticks = 5;
        visparams.tvscale = .06;
    end
    skipTV = ~isfield(res,'TVRes');
    
    if ~isfield(visparams,'tvscale')
        visparams.tvscale = .06;
    end
    BBcenter = [min(points)+ max(points)]/2;
    points = points - BBcenter;
    
    shape = polyshape(points);

    %% clean up data. Exclude TV nans and exclude past rlimlim
    % show MA data up to rlim
    if params.useNeckRlim && numel(res.ORes.keypoints.minRadNeck)~=0
        rlimused = res.ORes.keypoints.minRadNeck;
    else
        rlimused = params.rlimlim;
    end
    rlimArea = area(morphopen(shape, rlimused));
    keepindsMA = res.MARes.areas < rlimArea;
    % show TV data excluding nans or invalid solutions
    if ~skipTV; keepindsTV = ~isnan(res.TVRes.perims) & res.TVRes.perims < res.ORes.keypoints.perimMax*1.05; end;
    
    % compute left side tight region
    rt = max(res.MARes.props.maxneckrad, res.MARes.props.maxCirc2rad/2);
    inxy = mean(res.ORes.popens(2).Vertices);
    if isfield(res.ORes.keypoints, 'inxy')
        inxy = res.ORes.keypoints.inxy;
    end
    leftsidetightT = getArea1(shape, rt, inxy);
    % inx = mean(res.ORes.popens(2).Vertices);
    % [res.ORes.keypoints.inr res.MARes.props.maxneckrad] %0.3366    0.3304

    %% Visualize results
    fh=figure; ax0 = gca; hold all; set(gcf,'color','w');
    fh.Renderer = 'Painters';
    ylabel('Perimeter'); xlabel('Area'); title('Isoperimetric Profile Bound');
    if isfield(visparams,'name')
        title(visparams.name);
    end
    fh.OuterPosition = [1   1    700  350.2000];
    if isfield(visparams,'OuterPosition')
        fh.OuterPosition = visparams.OuterPosition;
    end
    bump = sqrt(var(res.ORes.perims))/2;
    xlim([res.ORes.keypoints.inrArea*.9 res.ORes.keypoints.areaMax*1.05]); 
    ylim([2*pi*res.ORes.keypoints.inr - 3*bump, res.ORes.keypoints.perimMax+2*bump]);
    if ~skipTV;
        plot([-10 0],[-10 0],'g','linewidth',2);
        plot([-10 0],[-10 0],'b','linewidth',2);
        plot([-10 0],[-10 0],'r','linewidth',2);
    else
        plot([-10 0],[-10 0],'g','linewidth',2);
        plot([-10 0],[-10 0],'b','linewidth',2);
    end
    if ~skipTV; plot(res.TVRes.areas(keepindsTV), res.TVRes.perims(keepindsTV), 'r.-','linewidth',2); end;
    plot(res.ORes.areas, res.ORes.perims, 'g-','linewidth',1.5); 
    plot(res.MARes.areas(keepindsMA), res.MARes.perims(keepindsMA), 'b-','linewidth',2)
    
    maxinrdir = 'right';
    if isfield(visparams,'maxinrdir')
        maxinrdir = visparams.maxinrdir;
    end
    xline(res.ORes.keypoints.inrArea,'k','linewidth',2,'label','Max inscribed circle','labelhorizontalalignment',maxinrdir)
    xline(res.ORes.keypoints.areaMax,'k','linewidth',2,'label','Max area','labelverticalalignment','bottom')
    if params.useNeckRlim
        xline(rlimArea,'k','linewidth',2,'label','Minimal Neck','labelverticalalignment','bottom','labelhorizontalalignment','left')
        xline(rlimArea,'m','linewidth',2)
    end
    if numel(res.MARes.props.isBigNeck)>0 && ~res.MARes.props.isBigNeck && exist('leftsidetightT','var') && numel(leftsidetightT)~=0
        xline(leftsidetightT,'k','linewidth',2,'label','Conservatively tight','labelhorizontalalignment','right')
        xline(leftsidetightT,'m','linewidth',2)
    end

    significantAreas = [linspace(res.ORes.keypoints.inrArea, rlimArea, visparams.nticks)];
    scale = .5;
    if isfield(visparams,'scale'); scale = visparams.scale; end;
    anisoscaleORat = fh.OuterPosition;
    anisoscaleO = anisoscaleORat(4)/anisoscaleORat(3);
    anisoscaleRat = get(gca,'DataAspectRatio');
    anisoscale = (anisoscaleRat(2)/anisoscaleRat(1)) / anisoscaleO;
    xlimb = xlim; ylimb = ylim;
    if params.useNeckRlim && (~isfield(visparams,'supressNeckShape') || ~visparams.supressNeckShape)
        iterlist = 1:numel(significantAreas);
    else
        iterlist = 1:numel(significantAreas)-1;
    end
    tvdownshift = 2;
    if isfield(visparams,'tvdownshift')
        tvdownshift = visparams.tvdownshift;
    end

    moupshift = 1;
    if isfield(visparams,'moupshift')
        moupshift = visparams.moupshift;
    end
    madownshift = 1;
    if isfield(visparams,'madownshift')
        madownshift = visparams.madownshift;
    end

    for i=iterlist
        areaI = significantAreas(i);

        if ~skipTV
            [~,ind2] = min(abs(res.TVRes.areas - areaI));
            perim2 = res.TVRes.perims(ind2);
            if ~isnan(perim2) % since mosek can sometimes fail, skip those indices
                ps2 = flipud(res.TVRes.popens(:,:,ind2));
                superimposeIndicatorOnFigure(fh, ps2, [areaI, perim2-tvdownshift*bump], visparams.tvscale);
            end
        end

        [~,ind0] = min(abs(res.ORes.areas - areaI));
        [~,ind1] = min(abs(res.MARes.areas - areaI));
        ps0 = res.ORes.popens(ind0);
        ps1 = res.MARes.popens(ind1);
        perim0 = res.ORes.perims(ind0);
        perim1 = res.MARes.perims(ind1);
        if exist('perim2','var') && isnan(perim2)
            perim2 = perim1;
        end
        if skipTV; perim2 = perim1; end;
        ps0T = transformShape(ps0, scale, [areaI, perim0+moupshift*bump], anisoscale);
        ps1T = transformShape(ps1, scale, [areaI, perim2-madownshift*bump], anisoscale);
        plot(ps0T,'facecolor','green','facealpha',.4);
        plot(ps1T,'facecolor','blue','facealpha',.8)
        
    end
    ps0T = transformShape(res.ORes.popens(end), scale, [res.ORes.keypoints.areaMax, res.ORes.keypoints.perimMax+moupshift*bump], anisoscale);
    ps1T = transformShape(res.MARes.popens(end), scale, [res.ORes.keypoints.areaMax, res.TVRes.perims(end-1)-madownshift*bump], anisoscale);
    plot(ps0T,'facecolor','green','facealpha',.4)
    plot(ps1T,'facecolor','blue','facealpha',.8)
    if ~skipTV
        ps2 = flipud(res.TVRes.popens(:,:,end));
        superimposeIndicatorOnFigure(fh, ps2, [res.ORes.keypoints.areaMax, res.TVRes.perims(end-1)-tvdownshift*bump], visparams.tvscale);
    end
    xlim(xlimb); ylim(ylimb);
    if ~isfield(visparams,'legend') || visparams.legend
        if ~skipTV
            lgd = legend({'Morph Open','Med Axis','TV'});
        else
            lgd = legend({'Morph Open','Med Axis'});
        end
    end
    fh.Children = fh.Children([end 1:end-1]);
    set(gca, 'Color', 'None')
    if ~isfield(visparams,'OuterPosition')
        fh.OuterPosition = [1   1    700  450];
    end
    
    if ~isfield(visparams,'legend') || visparams.legend
        if ~isfield(visparams,'lgdpos')
            lgd.Position = [0.2881    0.7420    0.1724    0.1681];
        else
            lgd.Position = visparams.lgdpos;
        end
    end
    
end