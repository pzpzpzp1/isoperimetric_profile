function superimposeIndicatorOnFigure(fh, ps2, center, scale)
    ps2 = indicator2redwhite(ps2);
    ax0 = fh.CurrentAxes;
    ax0pos = center;
    p0 = ax0.Position;
    xlim0 = ax0.XLim;
    ylim0 = ax0.YLim;
    xperc = (ax0pos(1)-xlim0(1))/(xlim0(2)-xlim0(1));
    yperc = (ax0pos(2)-ylim0(1))/(ylim0(2)-ylim0(1));
    perc = [xperc, yperc];
    axipos = [p0(1)+p0(3)*perc(1), p0(2)+p0(4)*perc(2)];
    
    axiposLL = axipos - [scale/2 scale/2];
    axi = axes('position',[axiposLL(1) axiposLL(2) scale scale],'visible','on');
    imshow(ps2);
    set(fh,'CurrentAxes',ax0);
end