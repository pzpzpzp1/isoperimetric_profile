clear all; close all;

files = dir('../*/cached*.mat');
envelopes = {};
names = {};
pointses = {};
for i=1:numel(files)
    load([files(i).folder '\' files(i).name]);
    try;
        envelopes{end+1} = extractAggregateIPEnvelope(points, res, .001);
        names{end+1} = files(i).name;
        pointses{end+1} = points;
    catch; 
        % might fail if TV part or O part wasn't saved or data is on older version of
        % schema. There's enough data regardless, so ignore.
    end;
end

f1 = figure; hold all; set(gcf,'color','w'); 
xlabel('Scaled Areas'); ylabel('Scaled Perimeters');
title('Aggregated Isoperimetric Profile Bounds');
xlim([0 1]);
pmax = -1;
colors = distinguishable_colors(numel(envelopes),'w');
for i=1:numel(envelopes)
    
    env = envelopes{i}; 
    name = names{i};

    % scale so all areas are 1. perimeters are sqrt scaled.
    area = max(env.Vertices(:,1));
    env.Vertices(:,1) = env.Vertices(:,1) * 1/area;
    env.Vertices(:,2) = env.Vertices(:,2) * sqrt(1/area);
    
    if max(env.Vertices(:,2))>pmax; 
        pmax = max(env.Vertices(:,2)); 
    end
    
    plt = plot(env,'edgecolor','none','linewidth',1,'facecolor',colors(i,:));
    plt.EdgeColor = plt.FaceColor;
    plt.FaceAlpha = 1;
end
ylim([0 pmax]);
legend(names);
f1.OuterPosition = [1   1    750  450];
exportgraphics(f1,'agg_IP_bound.pdf','contenttype','vector');


%% generate figures of all domains tested;
for i=1:numel(envelopes)
    name = names{i};
    
    ft = figure; axis equal; hold all; set(gcf,'color','w'); axis off;
    pts = pointses{i}; pts = pts - (min(pts)+max(pts))/2;
    pts = pts/max(max(pts));
    plot(polyshape(pts),'facecolor','k','facealpha',1);
    exportgraphics(ft,[names{i}(7:end-4) '_dom.pdf'],'contenttype','vector');
    close(ft);
end