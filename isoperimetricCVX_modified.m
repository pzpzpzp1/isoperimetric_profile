function [result,val] = isoperimetricCVX_modified(mask, percent, dxs)

if nargin < 3
    dxs = [1 1];
end

% A few special cases
if percent == 0
    result = zeros(size(mask));
    val = isotropicTotalVariation(result);
    return
end

if percent == 1
    result = mask;
    val = isotropicTotalVariation(result);
    return;
end

% Number of pixels to be used in the isoperimetric profile
pixelSum = percent * sum(mask(:));
mask = double(mask);
dx = dxs(1); dy = dxs(2);
cvx_begin
    cvx_solver mosek
    % cvx_precision best
    
    variable x(size(mask,1),size(mask,2))
  
    dfdx1 = (x(:,3:end)-x(:,1:end-2))/(2*dx);
    dfdx = (dfdx1(1:end-2,:) + dfdx1(2:end-1,:) + dfdx1(3:end,:))/3;
    dfdy1 = (x(3:end,:)-x(1:end-2,:))/(2*dy);
    dfdy = (dfdy1(:,1:end-2) + dfdy1(:,2:end-1) + dfdy1(:,3:end))/3;
    v = norms([dfdx(:) dfdy(:)],2,2);
    
    minimize sum(v)*dx*dy % perimeter
    subject to
        sum(x(:)) == pixelSum
        0 <= x <= mask
cvx_end

result = x;
val = cvx_optval;