close all; clear all;

districtfiles = dir('./ipvoronoi/data/All_states/AL_CON_Districts/PPoints/*.geojson');
statefiles = dir('./ipvoronoi/data/All_states/Points/*.geojson');
allfiles = [statefiles; districtfiles];

Ns = [5 10 20 40 100];
circleNs = [5 10 20 40 100];
names = {}; excepts = {};
for k=1:numel(allfiles)
    fname = [allfiles(k).folder '\' allfiles(k).name];
    points = loadGeojson(fname);
    for i=1:numel(Ns)
        for j=1:numel(circleNs)
            N = Ns(i);
            circleN = circleNs(j);
            display([i j k]);
            names{i,j,k}.N = N;
            names{i,j,k}.circleN = circleN;
            names{i,j,k}.fname = fname;
            try
                [fullperimeters{i,j,k}, fullareas{i,j,k}, xyrs{i,j,k}] = isoperim_profile_v2(points, 0, N, circleN);
            catch exception
                excepts(end+1).E = exception;
                excepts(end).fname = fname;
                excepts(end).N = N;
                excepts(end).circleN = circleN;
            end
        end
    end
end




