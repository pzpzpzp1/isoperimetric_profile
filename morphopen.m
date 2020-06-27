function [popen, perode] = morphopen(ps, rad)
    if numel(rad)==0 || rad==0 || numel(ps.Vertices)==0
        popen = ps;
        perode = ps;
        return;
    end
    
    perode = morpherode(ps,rad);
    popen = morphdilate(perode, rad);
    
end