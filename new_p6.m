function p6 = new_p6(tile)
    magfactor = 10;
    tile1 = imresize(tile, magfactor, 'bicubic');
    s = size(tile1, 1);

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
    tile0 = tile1(:, 1:x1);
    
    
    %right triangle inscribed into rectangle 
    %size(tile0) = [2y1 x x1]
    tile0 = tile0.*mask;
    
    %rotate tile1
    tile120 = imrotate(tile0, 120, 'bicubic');
    tile240 = imrotate(tile0, 240, 'bicubic');

    %trim the tiles manually, using trigonometric laws
    %AY NOTE: floor and round give us values that differ by 1 pix.
    %to trim right, we'll have to derive the trim value from 
    tile0 = [tile0 zeros(s, x1*2)];    
    delta = size(tile0, 2);
    
    %ideally we would've used  
    %delta = floor(sqrt(3)*s/2) OR round(sqrt(3)*s/2);   
    x120 = size(tile120, 2) - delta;
    y120 = size(tile120, 1) - y1;
    
    %size(tile120, 240) = [y1 x 3x1]
    tile120 = tile120(y120 + 1:end, x120 + 1:end);
    tile240 = tile240(1:y1, x120 + 1:end);
    
    %we have 3 tiles that will comprise
    %equilateral triangle together
    
    %glue them together by padding and smoothing edges (max instead of sum)
    %tile0 already padded
    tile120 = [zeros(y1, x1*3); tile120];
    tile240 = [tile240; zeros(y1, x1*3)];
    
    %size(tri) = [2y1 x 3x1]
    tri = max(max(tile0, tile120), tile240);
    %mirror_tri = fliplr(tri); --wrong! should be (fliplr(flipud(tri)))
    mirror_tri = rot90(tri, 2);
    
    %shifw w.slight overlap, 
    delta_pix = 3;
    shifted = [mirror_tri(y1 + 1 - delta_pix:end - delta_pix, :); mirror_tri(delta_pix + 1:y1 + delta_pix, :)];
    
    tile2 = max(tri, shifted);
    t2 = floor(0.5*size(tile2, 1));
    
    tile2_flipped = [tile2(t2 + 1:end, :); tile2(1:t2, :)]; 
    
    %size(tile3) = [2y1 x 6x1]
    tile3 = [tile2, tile2_flipped];
    p6 = imresize(tile3, 1/magfactor, 'bicubic');
end