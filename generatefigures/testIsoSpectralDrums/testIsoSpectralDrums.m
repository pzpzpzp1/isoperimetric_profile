clear all; close all;

D1 = [1 3;1 2;3 2;3 1;2 0;2 1;1 1;0 2; 1 3];
D2 = [0 3;1 3; 1 2; 2 2; 3 1;2 1;2 0;0 2; 0 3];
D1 = D1 - (min(D1)+max(D1))/2;
D2 = D2 - (min(D2)+max(D2))/2;

p1 = polyshape(D1);
p2 = polyshape(D2);
figure; 
subplot(1,2,1); hold all; axis equal; plot(p1); title('D1')
subplot(1,2,2); hold all; axis equal; plot(p2); title('D2')

[inr1, center1] = getMaxInscribedCircle(p1, 0);
[inr2, center2] = getMaxInscribedCircle(p2, 0);
pc1 = circlepoly_v3(inr1, center1, .01);
pc2 = circlepoly_v3(inr2, center2, .01);

%% plot medaxis 1
fig1 = visualizeMedAxis(D1);
plot(pc1,'linewidth',2,'facecolor','none','edgecolor','r');

fig2 = visualizeMedAxis(D2);
plot(pc2,'linewidth',2,'facecolor','none','edgecolor','r');