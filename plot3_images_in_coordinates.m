function plot3_images_in_coordinates(im, coordinates)
% TO plot images in a 3D plot, each one in a certain coordinate.
%
% INPUTS:
%       - im: image matrix of (N x M x CH),  of CH channels (it should be 1 or 3), 
%                     N by M pixels.
%       - coordinates: R x 3 matrix -> (x,y,z) coordinates for each img.
% _________________________________________________________________________

    
    %size of the images
    imr = imresize(im,.03);
    sizeim = size(imr);
    szx=20;
    szy=sizeim(1);
    szz=sizeim(2);
    szz=szy;
    minx=inf;
    miny=inf;
    minz=inf;
    maxx=-inf;
    maxy=-inf;
    maxz=-inf;
    

        hold on;      %# Add to the plot
        %xlabel(xlabelstr);
        %ylabel(ylabelstr);
        %zlabel(zlabelstr);
        img = im;  %# Load a sample image
        img(:,:,1)=flipud(img(:,:,1));
        img(:,:,2)=flipud(img(:,:,2));
        img(:,:,3)=flipud(img(:,:,3));

     
        xImage = [coordinates(1)-szx coordinates(1)+szx;coordinates(1)-szx coordinates(1)+szx];   %# The x data for the image corners
        yImage = [coordinates(2)-szy coordinates(2)+szy;coordinates(2)-szy coordinates(2)+szy];   %# The y data for the image corners
        zImage = [coordinates(3)-szz coordinates(3)-szz;coordinates(3)+szz coordinates(3)+szz];   %# The z data for the image corners
        surf(xImage,yImage,zImage,...    %# Plot the surface
             'CData',img,...
             'FaceColor','texturemap');
%         grid on
%         axis equal
 
        if minx>xImage(1,1)
            minx=xImage(1,1);
        end
 
        if miny>yImage(1,1)
            miny=yImage(1,1);
        end
        if minz>zImage(1,1)
            minz=zImage(1,1);
        end
        if maxx<xImage(2,2)
            maxx=xImage(2,2);
        end
        if maxy<xImage(2,2)
            maxy=xImage(2,2);
        end
        if maxz<zImage(2,2)
            maxz=zImage(2,2);
        end
    
       view([1 0 0]);

end