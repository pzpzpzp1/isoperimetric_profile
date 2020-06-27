function points = loadGeojson(fname)
    if nargin==0
        fname = [IPbasedir '\ipvoronoi\data\All_states\Points\05_t2_1.geojson'];
    end
    fid = fopen(fname);
    raw = fread(fid,inf);
    str = char(raw');
    fclose(fid);
    data = jsondecode(str);
    points = data.features.geometry.coordinates;
    if strcmp(class(points),'cell')
        points = points{end};
    end
    points = flipud(points);
    points = points - mean(points);
end