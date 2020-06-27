function [xys, radii] = getEdgeEdgePoints(radii_start_and_end, startpointc, endpointc, N)
startpoint = permute(startpointc(:,[1 2]),[2 3 1]);
endpoint = permute(endpointc(:,[1 2]),[2 3 1]);

%% sample points. Eventually could use arclength sampling, but for now, just uniform N samples in x will do.
t = linspace(0,1,N);
xys = startpoint.*(1-t) + endpoint.*(t);
radii = permute(radii_start_and_end(:,1).*(1-t) + radii_start_and_end(:,2).*(t),[3 2 1]);

end