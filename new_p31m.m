function p31m = new_p31m(tile)

    magfactor = 10;
    tile0 = imresize(tile, magfactor, 'nearest');
    s = size(tile0, 1);

    x1 = round(sqrt(3)*s/6);
    y1 = round(s/2);
    
    %vetrices of the triangle (closed polygon => four points)
    mask_x1 = [0 x1 0 0];
    mask_y1 = [0 y1 2*y1 0];

    %half of the mask
    %reflect and concatenate, to get the full mask:   
    
    mask_half = poly2mask(mask_x1, mask_y1, y1, x1);
    mask = cat(1, mask_half, flipud(mask_half));
    
    %copy and cut the tile 
    tile1 = tile0(:, 1:x1);
    
    %right triangle inscribed into rectangle 
    %size(tile1) = [2y1 x x1]
    tile1 = tile1.*mask;
    
    %%F
%     tile1 = place_letter(tile1);
    
    %rotate the tile
    tile120 = imrotate(tile1, 120);
    tile240 = imrotate(tile1, 240);
    
    %trim the tiles manually, using trigonometric laws
    %AY NOTE: floor and round give us values that differ by 1 pix.
    %to trim right, we'll have to derive the trim value from 
    tile1 = [tile1 zeros(s, x1*2)];    
    delta = size(tile1, 2);
    
    %ideally we would've used  
    %delta = floor(sqrt(3)*s/2) OR round(sqrt(3)*s/2);   
    x120 = size(tile120, 2) - delta;
    y120 = size(tile120, 1) - y1;

    %size(tile120, tile240) = [2y1 3x1]    
    tile120 = tile120(y120 + 1:end, x120 + 1:end);
    tile240 = tile240(1:y1, x120 + 1:end);
    
    %we have 3 tiles that will comprise
    %equilateral triangle together
    
    %glue them together by padding and smoothing edges (max instead of sum)
    %tile1 already padded
       
    tile120 = [zeros(y1, x1*3); tile120];
    tile240 = [tile240; zeros(y1, x1*3)];
    
    %size(tri) = [2y1 3x1]
    tri = max(max(tile1, tile120), tile240);
    mirror_tri = fliplr(tri);
    
    %use shift overlap, to smooth the edges
    delta_pix = 3;
    shifted = [mirror_tri(y1 + 1 - delta_pix:end - delta_pix, :); mirror_tri(delta_pix + 1:y1 + delta_pix, :)];
    tile2 = max(shifted, tri);

    %size(tile3) = [2y1 6x1]
    tile3 = [tile2, fliplr(tile2)];
    p31m = imresize(tile3, 1/magfactor, 'nearest');   
end
