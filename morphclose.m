function [pclose, pdilate] = morphclose(ps, rad)
    if numel(rad)==0 || rad==0 || numel(ps.Vertices)==0
        pdilate = ps;
        pclose = ps;
        return;
    end
    
    pdilate = morphdilate(ps, rad);
    pclose = morpherode(pdilate, rad);
end