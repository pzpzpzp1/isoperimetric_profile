function results = aggregateProcessing(points, params)

    if nargin==0
        points = loadCachedShape();
        
        params.dx = .1
        params.dc = .1
        params.rlimlim = mean(max(points) - min(points))/70;
        params.nHeightPixels = 10;
        params.Visualize = 1;
        params.useNeckRlim = 1;
        params.debug = 0;
        params.nMorphOpenSamples = 100;
    end
    if ~isfield(params,'postprocess')
        params.postprocess = 0;
        params.postprocessclosingradius = nan;
    end
    dx = params.dx;
    dc = params.dc;
    rlimlim = params.rlimlim;
    nHeightPixels = params.nHeightPixels;
    Visualize = params.Visualize;
    useNeckRlim = params.useNeckRlim;
    debug = params.debug;
    nMorphOpenSamples = params.nMorphOpenSamples;
    
    %% extract some basic info
    BBcenter = [min(points)+ max(points)]/2;
    points = points - BBcenter;
    shape = polyshape(points); 
    totalArea = area(shape); 
    totalPerimeter = perimeter(shape);
    if Visualize
        fh = visualizeMedAxis(points);
    end
    
    %% process morph open bound
    tic; [fullperimeters0, fullareas0, keypoints, popens0, perodes0] = morphOpenBound(points, Visualize, [nMorphOpenSamples, dx]);
    tmopen = toc;
    ORes.perims = fullperimeters0;
    ORes.areas = fullareas0;
    ORes.keypoints = keypoints;
    ORes.popens = popens0;
    ORes.perodes = perodes0;
    ORes.time = tmopen;
    results.ORes = ORes;
    
    %% process medial axis bound
    if useNeckRlim && numel(keypoints.minRadNeck)~=0
        rlim = max(keypoints.minRadNeck*.95, rlimlim); 
    else
        rlim = rlimlim; 
    end
    maparams.dc = dc;
    maparams.rlim = rlim;
    maparams.dx = dx;
    maparams.postprocess = params.postprocess;
    maparams.postprocessclosingradius = params.postprocessclosingradius;
    tic; [fullperimetersUpperBound2, fullareasMA2, popens2, props] = isoperim_profile_v3(points, Visualize, maparams, debug);
    tMA2 = toc;
    MARes.perims = fullperimetersUpperBound2;
    MARes.areas = fullareasMA2;
    MARes.popens = popens2;
    MARes.props = props;
    MARes.time = tMA2;
    results.MARes = MARes;
    
    %% process TV bound
    startAreaPercentage = keypoints.inrArea/totalArea;
    if keypoints.areaMinNeck/keypoints.areaMax ~= 1
        MidAreaSampleNumber = 15;
        EndAreaSampleNumber = 6;
        percentAreasMid = linspace(startAreaPercentage, keypoints.areaMinNeck/keypoints.areaMax, MidAreaSampleNumber);
        percentAreasEnd = linspace(keypoints.areaMinNeck/keypoints.areaMax, 1, EndAreaSampleNumber);
        percentAreas = [percentAreasMid percentAreasEnd(2:end)];
    else
        MidAreaSampleNumber = 20;
        percentAreas = linspace(startAreaPercentage, 1, MidAreaSampleNumber);
    end
    
    if ~isfield(params,'skipTV') || ~params.skipTV
        tic; [fullperimetersLowerBound1, fullareasTV1, images1] = getTVProfile(points, nHeightPixels, 'cvx', percentAreas);
        tTVcvx = toc;
        TVRes.perims = [0 fullperimetersLowerBound1'];
        TVRes.areas = [0 fullareasTV1];
        TVRes.popens(:,:,1) = images1(:,:,1)*0;
        TVRes.popens(:,:,2:size(images1,3)+1) = images1;
        TVRes.time = tTVcvx;
        results.TVRes = TVRes;
    end
end