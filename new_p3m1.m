function p3m1 = new_p3m1(tile)

    magfactor = 10;
    tile1= imresize(tile, magfactor, 'bicubic');
    height = size(tile1, 1);
    %% make equlateral triangular tile:
    
    %We'll need to create equlateral triangle mask with side length = s.
    %Create rectangle with sides s/2, s*sqrt(3)/2
    %according to Pythagorean theorem
    
    width = min(round(height*sqrt(3)), size(tile1, 2));

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
    tileR1 = imrotate(tile1_mirror, 240, 'bicubic');
    %AY: I directly cut the tiles, because trim will
    %return slightly different size
    
    t_r1x = size(tileR1, 1);
    tileR1 = tileR1(t_r1x - height + 1:end, 1:width);
    
    %AY: rotating mirrored tile(as opposed to tileR1) will cause less
    %border effects when we'll add it to two other tiles.
    tileR2 = imrotate(tile1_mirror, 120, 'bicubic');
    tileR2 = tileR2(1:height, 1:width);
    %%Assembling the tiles
    
    %We have 3 tiles with the same triangle rotated 0, 120 and 240
    %pad them and put them together
    zero_tile = zeros(y1, width);
    tile2 = [zero_tile; tile0; zero_tile];
    tileR1 = [zero_tile; zero_tile; tileR1];
    tileR2 = [tileR2; zero_tile; zero_tile];
    
    %Using max() will give us smoother edges, as opppose to sum()
    half1 = max(tile2, tileR1);
    half = max(half1, tileR2);
    
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
