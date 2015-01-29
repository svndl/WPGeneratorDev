function p6m = new_p6m(tile)

    magfactor = 10;
    tile1 = imresize(tile, magfactor, 'bicubic');
    height = size(tile1, 1);
    
    width = round(height/sqrt(3));

    %% fundamental region is right triangle with angles (30, 60)
    
    %vetrices of the triangle (closed polygon => four points)
    mask_x1 = [0 width 0 0];
    mask_y1 = [0 height height 0];
    
    %half of the mask
    %reflect and concatenate, to get the full mask:   
    mask = poly2mask(mask_x1, mask_y1, height, width);
    
    %right triangle inscribed into rectangle 
    tile0 = tile1(:, 1:width).*mask;
    
    %size(tile0) = [height x width]   
    tile0 = [tile0; flipud(tile0)]; 
    
    %rotate tile1
    tile120 = imrotate(tile0, 120, 'bilinear');
    tile240 = imrotate(tile0, 240, 'bilinear');

    %trim the tiles manually, using trigonometric laws
    %NOTE: floor and round give us values that differ by 1 pix.
    %to trim right, we'll have to derive the trim value from 
    tile0 = [tile0 zeros(2*height, width*2)];    
    delta = size(tile0, 2);
    
    %ideally we would've used  
    %delta = floor(sqrt(3)*s/2) OR round(sqrt(3)*s/2);   
    x120 = size(tile120, 2) - delta;
    y120 = size(tile120, 1) - height;
    
    %size(tile120, 240) = [y1 x 3x1]
    tile120 = tile120(y120 + 1:end, x120 + 1:end);
    tile240 = tile240(1:height, x120 + 1:end);
    
    %we have 3 tiles that will comprise
    %equilateral triangle together
    
    %glue them together by padding and smoothing edges (max instead of sum)
    %tile0 already padded
    tile120 = [zeros(height, width*3); tile120];
    tile240 = [tile240; zeros(height, width*3)];
    
    %size(tri) = [2y1 x 3x1]
    tri = max(max(tile0, tile120), tile240);
    mirror_tri = rot90(tri, 2);
    
    %shifw w.slight overlap, 
    delta_pix = 3;
    shifted = [mirror_tri(height + 1 - delta_pix:end - delta_pix, :); mirror_tri(delta_pix + 1:height + delta_pix, :)];
    
    tile2 = max(tri, shifted);
    t2 = floor(0.5*size(tile2, 1));
    
    %tile2_flipped = [tile2(t2 + 1:end, :); tile2(1:t2, :)];
    tile2_flipped = rot90(tile2, 2); 
    
    %size(tile3) = [2y1 x 6x1]
    tile3 = [tile2, tile2_flipped];
    p6m = imresize(tile3, 1/magfactor, 'bicubic');
end
