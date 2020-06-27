
[img, map, transp] = imread('mz-01.png');
points = userinput_points('cachedshape.lines',transp); 
save('mozpoints.mat','points')