function T = extractRuntimeInfo(res, points, params, name)
    numVerts = size(points,1);

    dc = params.dc;
    rlim = params.rlimlim;
    rn = res.ORes.keypoints.minRadNeck;
    rm = res.ORes.keypoints.maxRadNeck;
    MAtime = res.MARes.time;

    if numel(rn)==0
        rn = inf;
        rm = inf;
    end

    nOSamples = params.nMorphOpenSamples;
    Otime = res.ORes.time;

    skipTV = isfield(params,'skipTV') && params.skipTV;
    if ~skipTV
        nTVSamples = numel(res.TVRes.perims);
        npix = params.nHeightPixels;
        TVTime = res.TVRes.time;
    else
        nTVSamples = 0;
        npix = 0;
        TVTime = 0;
    end
    T = table(numVerts, dc, rlim, rn, rm, MAtime, nOSamples, Otime, nTVSamples, npix, TVTime, {name});    
    T.Properties.VariableNames = {'Number of Vertices','dc','rlim','rn','rm','MA time', 'Number of Open Samples', 'Open time', 'Number of TV Samples', 'Vertical Pixels Resolution','TV time','Name'}; 
end









