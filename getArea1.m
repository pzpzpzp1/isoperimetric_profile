function z = getArea1(ps, radius, inxy)
    if numel(radius)==0
        z = [];
        return;
    end
    
    rplus = radius+.0005;
    [popen, perode] = morphopen(ps, rplus);
    
    perodeRegions = regions(perode);
    isIn = zeros(numel(perodeRegions),1) == 1;
    for i=1:numel(perodeRegions)
        isIn(i) = isinterior(perodeRegions(i),inxy(1),inxy(2));
    end
    assert(any(isIn)); % if this is false, then check if you centered the shape based on bounding box center first in all computations.
    subR1 = morphdilate(perodeRegions(isIn), rplus);
    z = area(subR1);
    %{
    figure; hold all; axis equal; 
    plot(ps);
    plot(popen);
    plot(subR1);
    plot(perode);
    plot(inxy(1),inxy(2),'r.')
    %}
end










