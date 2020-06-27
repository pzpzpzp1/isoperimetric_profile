function points = scaleAndCenter(points)
    BBcenter = (min(points)+max(points))/2;
    points = points - BBcenter;
    points = points / max(max(points));
end
