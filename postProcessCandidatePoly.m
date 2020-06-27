% given a candidate inscribed poly of some area, it can use some cleanup
% due to little polygon artifacts from using all these polyshape unions.
% 1. no holes in the candidate. There's never a reason to have holes (unless the domain has holes. that'll be handled later)
% 2. closed for some radius. all free boundaries are supposed to be circular so they don't change under closing. Closing can increase area outside boundary, but we can handle that later.
% 3. ok now it's later, just intersect the candidate with the original poly: base.
% One caveat: If the original domain has boundary components very close to each other, and closing radius is too big so that candidate on one side bleeds into candidate on the other side. solution, choose a small closing radius :/ sadly that's another parameter.
function poly = postProcessCandidatePoly(base, poly, closeRad)
    poly = rmholes(poly);
    poly = morphclose(poly, closeRad);
    poly = intersect(poly, base);
end