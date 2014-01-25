function p3 = new_p3(tile)

   magfactor = 5;
   %resulting tile size will be
   s0 = size(tile, 1);
%    tile1 = imresize(tile, [3*s0*magfactor/2 round(s0*magfactor*sqrt(3)/2)], 'nearest');

   newSize = 3*s0*magfactor/2;
   newWidth = round(s0*magfactor*sqrt(3)/2);
   tile1 = imresize(tile, [newSize newSize], 'nearest');   
   tile1 = tile1(:, 1:newWidth);
   
   %AY Jan 7, 14: NEW Code. Instead of tilting the tile, use a mask to cut the rhombus    
   %make a mask
   
   height = size(tile1, 1);
   width = size(tile1, 2);
   s1 = round(width/sqrt(3));
   x = [0 width width 0 0];
   y = [0 s1 height 2*s1 0];
   mask = poly2mask(x, y, height, width);
   tile0 = mask.*tile1;
   s = 2*s1;
    
   %AY Jan 7, 14: OLD Code    
    %tilt the tile into rhombus, inscribe the rhombus into rectangle
%    tile1 = imresize(tile, magfactor, 'bicubic');    
%    s = size(tile1, 1);   
%    length = s*1.5;
%    width = round(sqrt(3)*s/2);
%    tile0 = zeros(length, width);
%    for i = 1:width 
%        tile0(1 + floor(i/sqrt(3)):s + floor(i/sqrt(3)), i) = tile1(:, i);
%    end


    %rotate rectangle by 120, 240 degs
   tile120 = imrotate(tile0, 120, 'bilinear');
   tile240 = imrotate(tile0, 240, 'bilinear');
   
   %%manually trim the tiles:
   
   %tile120 should have the same size as rectangle: [heigh x width]
   tile120 = tile120(1:height, (floor(0.5*width) + 1):(floor(0.5*width) + width));
   
   %tile240 should have the size [s x 2*width]
   %find how much we need to cut from both sides
   diff = round(0.5*(size(tile240, 2) - 2*width));
   tile240 = tile240((round(0.25*s) + 1):(round(0.25*s) + s), (diff + 1):(2*width + diff));
   
   %%Start to pad tiles and glue them together
   %Resulting tile will have the size [2*length x 2* width]
   
   two_thirds1 = [tile0 tile120];
   two_thirds2 = [tile120 tile0];
   
   two_thirds = [two_thirds1; two_thirds2];
 
   
   %lower half of tile240 on the top, zero-padded to [length x 2 width]
   one_third11 = [tile240(0.5*s + 1:end, :); zeros(s, 2*width)];
   
   %upper half of tile240 on the bottom, zero-padded to [length x 2 width]
   one_third12 = [zeros(s, 2*width); tile240(1:0.5*s, :)];
    
   %right half of tile240 in the middle, zero-padded to [length x 2 width]
   one_third21 = [zeros(s, width); tile240(:, width + 1:end); zeros(s, width)];
   
   %left half of tile240in the middle, zero-padded to [length x 2 width]
   one_third22 =  [zeros(s, width); tile240(:, 1:width); zeros(s, width)];
   
   %cat them together
   one_third1 = [one_third11; one_third12];
   one_third2 = [one_third21 one_third22];

   %%glue everything together, shrink and replicate 
   one_third = max(one_third1, one_third2);
   
   %size(whole) = [2xlength 2xwidth]
   whole = max(two_thirds, one_third);
   p3 = imresize(whole, 1/magfactor, 'nearest');  
end