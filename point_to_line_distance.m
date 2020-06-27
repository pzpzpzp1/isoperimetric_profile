function d = point_to_line_distance(pt, v1, v2)

pt = [pt zeros(size(pt,1),1)];
v1 = [v1 zeros(size(v1,1),1)];
v2 = [v2 zeros(size(v2,1),1)];

a = v1 - v2;
b = pt - v2;
d = sqrt(sum(cross(a,b,2).^2,2)) ./ sqrt(sum(a.^2,2));