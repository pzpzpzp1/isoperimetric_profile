% converts 1 indicator to red.
% converts 0 indicator to white.
% rest is pink.
function z = indicator2redwhite(ps2)
    z(:,:,1) = 0*ps2+1;
    z(:,:,2) = 1-ps2;
    z(:,:,3) = 1-ps2;
end