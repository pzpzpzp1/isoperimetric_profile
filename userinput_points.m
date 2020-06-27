function xys=userinput_points(fname, image)
    
    fh=figure; hold all; axis equal;
    title('click to add points, middle click to stop.'); 
    
    if nargin >= 2
        image = flipud(image);
        imagesc(image);
        xlim([1 size(image,2)]); ylim([1 size(image,1)]); 
    else
        xlim([-5 5]); ylim([-5 5]); 
    end

    set(gca,'xtick',[]);
    set(gca,'ytick',[]);
    % axis off;
    c=1;
    xys = [];
    while c~=2
        [x,y,c] = ginput(1);
        if c==1
            xys = [xys; x y];
            try; delete(gc); catch; end;
            gc = plot(xys(:,1),xys(:,2),'r.-','linewidth',2,'markersize',20);
            if nargin >=2; xlim([1 size(image,2)]); ylim([1 size(image,1)]); end;
            drawnow;
        end
    end
    xys = [xys; xys(1,:)];
    
    pause(.001);
    close(fh);
    
    if size(xys,1)<=3
        error('Not enough vertices given!')
    end
    
    if nargin~=0
        fid = fopen(fname,'w');
        fprintf(fid,'%f %f \n',xys');
        fclose(fid);
    end
end