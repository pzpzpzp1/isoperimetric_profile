% pulls cached data from the repository. if no runtimes are shown, then you have to run the figure generating code to generate data first. Though runtimes will be relative to the computer that ran it.
files = dir('../*/cached*.mat');

T = table;
for i=1:numel(files)
    load([files(i).folder '\' files(i).name])
    T = [T;extractRuntimeInfo(res, points, params, files(i).name)];
end

% FOR FULL RUNTIME PAUSE RUNTIMEs HERE.


% pick a good coverage of runtimes to make table and reformat for easy
% pasting
T2 = T(:,[6 end 8 11 1 9 10 4]);
[~,reorder] = sort(T2.("Number of Vertices"));
T2 = T2(reorder,:);
selectedNames = {'cachedStar.mat', 'cachedJ.mat', 'cachedJellyfish.mat', 'cachedC.mat', 'cachedF.mat','cachedCandy.mat'};
inds = ismember(T2.Name, selectedNames);
finalT = T2(inds,[2 5 1 3 4 end]);

