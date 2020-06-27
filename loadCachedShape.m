function points = loadCachedShape(fname)
    if nargin == 0
        fname = 'cachedshape.lines';
    end
    fid = fopen(fname,'r'); 
    if fid == -1
        error('No cached shape available!');
    end
    points = fscanf(fid,'%f'); 
    fclose(fid); 
    points = reshape(points,2,[])'; 
end