function p3m1 = new_p3m1(tile)

    magfactor = 10;
    tile1= imresize(tile, magfactor, 'bicubic');
    height = size(tile1, 1);
    
    
    %% fundamental region is equlateral triangle with side length = height 
    width = round(0.5*height*sqrt(3));

    y1 = round(height/2);
   
    %vetrices of the triangle (closed polygon => four points)
    mask_x1 = [0  width 0 0];
    mask_y1 = [0 y1 y1 0];
    
    %Create half of the mask
    %reflect and concatenate, to get the full mask:   

    mask_half = poly2mask(mask_x1, mask_y1, y1, width);
    mask = cat(1, mask_half, flipud(mask_half));
    
    %equilateral triangle inscribed into rectangle 
    tile0 = tile1(:, 1:width).*mask;
 
    %%continue to modify the tile
    
    %reflect and rotate
    tile1_mirror = fliplr(tile0);
    tile240 = imrotate(tile1_mirror, 240, 'bilinear');
    %AY: I directly cut the tiles, because trim will
    %return slightly different size
    
    t_r1x = size(tile240, 1);
    tile240 = tile240(t_r1x - height + 1:end, 1:width);
    
    %AY: rotating mirrored tile(as opposed to tileR1) will cause less
    %border effects when we'll add it to two other tiles.
    tile120 = imrotate(tile1_mirror, 120, 'bilinear');
    tile120 = tile120(1:height, 1:width);
    %%Assembling the tiles
    
    %We have 3 tiles with the same triangle rotated 0, 120 and 240
    %pad them and put them together
    zero_tile = zeros(y1, width);
    tile2 = [zero_tile; tile0; zero_tile];
    tile240 = [zero_tile; zero_tile; tile240];
    tile120 = [tile120; zero_tile; zero_tile];
    
    %Using max() will give us smoother edges, as opppose to sum()
    half1 = max(tile2, tile240);
    half = max(half1, tile120);
    
    %By construction, size(whole) = [4*y1 2*x1]
    whole = [half, fliplr(half)];
    
    %Shifting by 2 pix (as oppose to zero), we'll smoothly glue tile together. 
    %Set delta_pix value to zero and see the difference
    delta_pix = 2;
    topbit = [whole((3*y1 + 1 - delta_pix):4*y1 -delta_pix, (width + 1):2*width) whole((3*y1 + 1- delta_pix):4*y1 - delta_pix, 1:width)];
    botbit = [whole(1 + delta_pix:y1 + delta_pix, (width + 1):2*width) whole(1 + delta_pix:y1 + delta_pix, 1:width)];
    
    whole(1:y1, :) = max(whole(1 + delta_pix:y1 + delta_pix, :), topbit);
    whole(3*y1 + 1:4*y1, :) = max(whole(3*y1 + 1 - delta_pix:4*y1 - delta_pix, :), botbit);
            
    %cutting middle piece of tile
    mid_tile = whole(y1 + 1:3*y1, 1:width);
    %reflecting middle piece and glueing both pieces to the bottom
    %size(bigTile)  = [6*y1 2*x1]    
    bigTile = [whole;  fliplr(mid_tile), mid_tile];
    p3m1 = imresize(bigTile, 1/magfactor, 'bicubic');
end
