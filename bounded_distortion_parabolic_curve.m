f = @(x)x.^2/10;
dx = .001;
t = 0:dx:1000;
x=sqrt(t); % edgelens converge. max edgeln -> sqrt(dx). min(edgelen) -> dx
% x=(t); % edgelens diverge
% x=log(t); % edglens go to 0
y=f(x);

xy=[x;y];
elens = vecnorm(xy(:,2:end)-xy(:,1:end-1),2,1);
plot(elens)

plot(x,y,'.')
