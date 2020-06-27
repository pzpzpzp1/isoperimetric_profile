clear all; close all;

img = imread('star.png');
points = userinput_points('cachedshape.lines', img);
points = scaleAndCenter(points);
save('starPoints.mat','points');





