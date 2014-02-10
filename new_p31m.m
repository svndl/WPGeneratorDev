function p31m = new_p31m(tile)

    magfactor = 10;
    tile0 = imresize(tile, magfactor, 'bicubic');
    
    height = size(tile0, 1);
    width = round(0.5*height/sqrt(3));
    
    y1 = round(height/2);
    
    %% fundamental region is an isosceles triangle with angles(30, 120, 30)
    
    %vetrices of the triangle (closed polygon => four points)
    mask_x1 = [0 width 0 0];
    mask_y1 = [0 y1 height 0];

    %make half of the mask
    %reflect and concatenate, to get the full mask:   
    
    mask_half = poly2mask(mask_x1, mask_y1, y1, width);
    mask = cat(1, mask_half, flipud(mask_half));

    
    %size(tile0) = [height  width]
    tile0 = mask.*tile0(:, 1:width);
    
    %rotate the tile
    tile120 = imrotate(tile0, 120, 'bilinear');
    tile240 = imrotate(tile0, 240, 'bilinear');
    
    %trim the tiles manually, using trigonometric laws
    %NOTE: floor and round give us values that differ by 1 pix.
    %to trim right, we'll have to derive the trim value from 
    tile0 = [tile0 zeros(height, width*2)];    
    delta = size(tile0, 2);
    
    %ideally we would've used  
    %delta = floor(sqrt(3)*s/2) OR round(sqrt(3)*s/2);   
    x120 = size(tile120, 2) - delta;
    y120 = size(tile120, 1) - y1;

    %size(tile120, tile240) = [height 3width]    
    tile120 = tile120(y120 + 1:end, x120 + 1:end);
    tile240 = tile240(1:y1, x120 + 1:end);
    
    %we have 3 tiles that will comprise
    %equilateral triangle together
    
    %glue them together by padding and smoothing edges (max instead of sum)
    %tile1 already padded
       
    tile120 = [zeros(y1, width*3); tile120];
    tile240 = [tile240; zeros(y1, width*3)];
    
    %size(tri) = [height 3width]
    tri = max(max(tile0, tile120), tile240);
    mirror_tri = fliplr(tri);
    
    %use shift overlap, to smooth the edges
    delta_pix = 3;
    shifted = [mirror_tri(y1 + 1 - delta_pix:end - delta_pix, :); mirror_tri(delta_pix + 1:y1 + delta_pix, :)];
    tile2 = max(shifted, tri);

    %size(tile3) = [height 6width]
    tile3 = [tile2, fliplr(tile2)];
    p31m = imresize(tile3, 1/magfactor, 'bicubic');   
end
