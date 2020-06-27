clear all; close all;

files = dir('../*/cached*.mat');
names = {};
pointses = {};
MAEs = {};
MOEs = {};
TVEs = {};
for i=1:numel(files)
    load([files(i).folder '\' files(i).name]);
    try;
        names{end+1} = files(i).name;
        pointses{end+1} = points;
        MAEs{end+1} = res.MARes.popens;
        MOEs{end+1} = res.ORes.popens;
        TVEs{end+1} = res.TVRes.popens;
    catch; 
        % might fail if TV part or O part wasn't saved or data is on older version of
        % schema. There's enough data regardless, so ignore.
    end;
end

colors = spring;
for i=1:numel(names)
    fh = figure; axis equal; hold all; axis off; set(gcf,'color','w');
    shf = MOEs{i}(end);
    for j=1:numel(MOEs{i})
        shi = MOEs{i}(numel(MOEs{i}) - j + 1);
        percentage = area(shi)/area(shf);
        colind = ceil(percentage*size(colors,1)); if colind==0; colind = 1; end;
        plot(MOEs{i}(numel(MOEs{i}) - j + 1),'facecolor',colors(colind,:),'edgecolor','none','linewidth',.5,'facealpha',1);
        drawnow;
    end
    exportgraphics(fh,[names{i}(7:end-4) '_MO.png'],'contenttype','image','Resolution',500);
    close(fh);
    
    fh = figure; axis equal; hold all; axis off; set(gcf,'color','w');
    for j=1:numel(MAEs{i})
        shi = MAEs{i}(numel(MAEs{i}) - j + 1);
        percentage = area(shi)/area(shf);
        colind = ceil(percentage*size(colors,1)); if colind==0; colind = 1; end;
        plot(MAEs{i}(numel(MAEs{i}) - j + 1),'facecolor',colors(colind,:),'edgecolor','none','linewidth',.5,'facealpha',1);
        drawnow;
    end
    exportgraphics(fh,[names{i}(7:end-4) '_MA.png'],'contenttype','image','Resolution',500);
    close(fh);
end
