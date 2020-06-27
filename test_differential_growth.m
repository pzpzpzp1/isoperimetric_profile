%% test differential growth edge-edge type
clear; close all;
D = rand*3;
r1 = rand; r2 = rand;
theta = (linspace(0,2*pi,1000)'); theta = theta(1:end-1); cs = flipud([cos(theta) sin(theta)]);
ts=linspace(0,.001,30);
phi = acos((r1-r2)/D);
A0 = r1^2*phi;
P0 = 2*r1*phi;
A1 = r1^2*(pi-phi);
P1 = 2*r1*(pi-phi);
dt = D*sin(phi);

perims = zeros(numel(ts),1);
areas = zeros(numel(ts),1);
r = zeros(numel(ts),1);
for i=1:numel(ts)
    t=ts(i);
    r(i) = r2*t+(1-t)*r1;
    d(i) = dt*t;
    A(i) = A1 + d(i)*(r1+r(i)) + r(i)^2*phi;
    P(i) = P1 + 2*d(i) + 2*r(i)*phi;
    
    s(i) = convhull(union(polyshape(r1*cs), polyshape(r(i)*cs+[D*t 0])));
    areas(i) = area(s(i));
    perims(i) = perimeter(s(i));
end

figure; 
subplot(2,2,1); hold all; xlabel('areas'); ylabel('perimeters');
plot(areas, perims, 'g'); 
plot(A, P, 'r');
subplot(2,2,3); hold all; xlabel('areas'); ylabel('growth-rate');
slopes = (perims(2:end)-perims(1:end-1))./(areas(2:end)-areas(1:end-1));
title([num2str(slopes(1)) ' vs ' num2str(1/r1)])
plot(areas(1:end-1),slopes,'r')
subplot(1,2,2); hold all; axis equal; 
plot(s(1)); 
plot(s(floor(numel(ts)/2))); 
plot(s(end));
plot(polyshape(r2*cs + [D 0]))
title(['nochange: ' num2str(D+r2 < r1)])



%% test differential growth vert-vert type.
% growth must be to the right. 
clear; close all;
x0 = 3*randn;
H = rand*7;
ts=linspace(x0,x0+.0001,30);
theta = (linspace(0,2*pi,1000)'); theta = theta(1:end-1); cs = flipud([cos(theta) sin(theta)]);
assert(ts(end) > x0); 
r0 = sqrt(H^2+x0^2);
A0 = pi*r0^2;
P0 = 2*pi*r0;
phiL = atan2(H,abs(x0)); if x0<0; phiL = pi-phiL; end;

perims = zeros(numel(ts),1);
areas = zeros(numel(ts),1);
for i=1:numel(ts)
    t=ts(i);
    r(i) = sqrt(t^2+H^2);
    phiR = atan2(H,abs(t)); if t>0; phiR = pi-phiR; end;
    A(i) = pi*r0^2*(phiL/pi) + pi*r(i)^2*(phiR/pi) + H*(t-x0);
    P(i) = 2*pi*r0*(phiL/pi) + 2*pi*r(i)*(phiR/pi);
    
    s(i) = union(polyshape(r0*cs+[x0 0]), polyshape(r(i)*cs + [t 0]));
    areas(i) = area(s(i));
    perims(i) = perimeter(s(i));
end

figure; 
subplot(2,2,1); hold all; xlabel('areas'); ylabel('perimeters');
plot(areas, perims, 'g'); 
plot(A, P, 'r.');
subplot(2,2,3); hold all; xlabel('areas'); ylabel('growth-rate');
slopes = (perims(2:end)-perims(1:end-1))./(areas(2:end)-areas(1:end-1));
title([num2str(slopes(1)) ' vs ' num2str(1/r0)])
plot(areas(1:end-1),slopes,'r')
subplot(1,2,2); hold all; axis equal; 
plot(s(1)); plot(s(end));
plot(0,H,'r.');plot(0,-H,'r.');plot(x0,0,'k.');plot(ts(end),0,'b.')
title(['err ' num2str(norm(areas-A') + norm(perims-P'))])



%% test differential growth edge-vert type.
% grow in positive direction from x0 in real.
% numerical evidence says growth is always 1/r0 
clear; close all;
H = rand*3;
x0 = rand*20-10;
% H = 1.4;
% x0 = 1.2;
a = 1/(2*H); 
y=@(x)a*x.^2+H/2;
ts=linspace(x0,x0+abs(x0),100);
theta = (linspace(0,2*pi,2000)'); theta = theta(1:end-1); cs = flipud([cos(theta) sin(theta)]);
r0 = y(x0);
A0 = pi*r0^2;
P0 = 2*pi*r0;
perims = zeros(numel(ts),1);
areas = zeros(numel(ts),1);
figure; axis equal; hold all; yline(0); plot(0,H,'k.'); plot(x0,0,'r.'); plot(ts(end),0,'r.')
for i=1:numel(ts)
    t=ts(i);
    r(i) = y(t);
    
    plot(t,y(t),'k.')
    plot(x0,y(x0),'k.')
    %coeffs = polyfit([x0,t],[y(x0),y(t)],1); xL = x0-r0;    xR = t+r(i);    yL = coeffs(1)*xL + coeffs(2);    yR = coeffs(1)*xR + coeffs(2);
    %plot(xL,yL,'k.');    plot(xR,yR,'k.')
    
    m = (y(t)-y(x0))/(t-x0);
    phi = pi/2-(atan2((y(t)-y(x0)), (t-x0)));
    A1 = pi*r0^2*phi/(2*pi);
    A2 = pi*r(i)^2*(pi-phi)/(2*pi);
    A3 = (t-x0)*(y(t)+y(x0))/2;
    p0 = [0,H];
    p1 = [x0,y(x0)];
    p2 = [t,y(t)];
    %p012 = polyshape([p0;p1;p2]);
    e01 = p0-p1;
    e02 = p0-p2;
    e12 = p2-p1;
    e21 = p1-p2;
    phiL=acos(dot(e01,e12)/norm(e01)/norm(e12));
    phiR=acos(dot(e02,e21)/norm(e02)/norm(e21));
    A4 = pi*r0^2*(pi-phiL)/(2*pi);
    A6 = pi*r(i)^2*(pi-phiR)/(2*pi);
    A5 = norm(cross([e01 0],[e12 0]))/2;
    P1 = 2*pi*r0*(phi+(pi-phiL))/(2*pi);
    P2 = 2*pi*r(i)*((pi-phi)+(pi-phiR))/(2*pi);
    P3 = t-x0;
    A(i) = A1 + A2 + A3 + A4 + A5 + A6;
    P(i) = P1 + P2 + P3;
    
    tri = polyshape([x0,0; t,0; 0,H]);
    s(i) = union(union(polyshape(r0*cs+[x0 y(x0)]), polyshape(r(i)*cs + [t y(t)])),tri);
    %s123 = intersect(s(i), polyshape([xL,yL;xR,yR;xR,0;xL 0]));
    %plot(s123);
    plot(s(i));
    plot(tri);
    areas(i) = area(s(i));
    perims(i) = perimeter(s(i));
end

figure; 
subplot(2,2,1); hold all; xlabel('areas'); ylabel('perimeters');
plot(areas, perims, 'g'); 
plot(A, P, 'r.');
subplot(2,2,3); hold all; xlabel('areas'); ylabel('growth-rate');
slopes = (perims(3:end)-perims(1:end-2))./(areas(3:end)-areas(1:end-2));
title([num2str(slopes(1)) ' vs ' num2str(1/r0)])
plot(areas(2:end-1),slopes,'r.')
plot(areas,1./r,'g')
subplot(1,2,2); hold all; axis equal; 
plot(s(1)); plot(s(end));
yline(0); plot(0,H,'k.'); plot(x0,0,'r.'); plot(ts(end),0,'r.')
title(['err ' num2str(norm(areas(2:end)-A(2:end)') + norm(perims(2:end)-P(2:end)'))])

%% not strictly necessary but test differential growth of a single circle out of nothing
clear; close all;
rs = linspace(0,1,100);
theta = (linspace(0,2*pi,2000)'); theta = theta(1:end-1); cs = flipud([cos(theta) sin(theta)]);
areas = zeros(numel(rs),1);
perims = zeros(numel(rs),1);
for i=1:numel(rs)
    s(i) = polyshape(rs(i)*cs);
    areas(i) = area(s(i));
    perims(i) = perimeter(s(i));
    A(i)=pi*rs(i)^2;
    P(i)=2*pi*rs(i);
end

figure; 
subplot(2,2,1); hold all; xlabel('areas'); ylabel('perimeters');
plot(areas, perims, 'g'); 
plot(A, P, 'r.');
subplot(2,2,3); hold all; xlabel('areas'); ylabel('growth-rate');
slopes = (perims(3:end)-perims(1:end-2))./(areas(3:end)-areas(1:end-2));
title([num2str(slopes(1)) ' vs ' num2str(inf)])
plot(areas(2:end-1),slopes,'r.');
plot(areas,1./rs,'g');
subplot(1,2,2); hold all; axis equal; 
plot(s(2)); plot(s(end));
title(['err ' num2str(norm(areas(2:end)-A(2:end)') + norm(perims(2:end)-P(2:end)'))])














