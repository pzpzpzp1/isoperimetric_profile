function fh = animate(popens1)
    assert(strcmp(class(popens1),'polyshape'));
    if numel(popens1)==0
        return;
    end
    
    fh = figure; 
    axis equal;
    hold all;
    drawnow;
    plot(popens1(end),'facecolor','blue');
    for i=1:numel(popens1)
        try; delete(ch); catch; end;
        ch = plot(popens1(i),'facecolor','green'); 
        drawnow; 
    end
end