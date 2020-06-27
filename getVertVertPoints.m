function [xys, radii] = getVertVertPoints(segment_point_cells, startpointc, endpointc, N)
startpoint = permute(startpointc(:,[1 2]),[2 3 1]);
endpoint = permute(endpointc(:,[1 2]),[2 3 1]);
V1 = permute(reshape(cell2mat(segment_point_cells(:,1)),2,[])',[2 3 1]);

%% sample points. Eventually should use arclength sampling, but for now, just uniform N samples in x will do.
t = linspace(0,1,N);
xys = startpoint.*(1-t) + endpoint.*(t);
radii = vecnorm(xys-V1,2,1);

end