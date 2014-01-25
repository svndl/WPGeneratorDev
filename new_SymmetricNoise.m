function image = new_SymmetricNoise(type, N, n, optTexture)

%Generates single wallaper group image
%--wptype defines the wallpaper group
%--N is the size of the output image. For complex groups
%returned image will have size at least NxN.
%--n the size of repeating pattern for all groups.
    
    %default 
    if(nargin < 4) 
        texture = filterTile(rand(n));
    else
        minDim = min(size(optTexture));
        %stretch user-defined texture, if it is too small for sampling
        if(minDim < n)
            ratio = round(n/minDim);
            optTexture = imresize(optTexture, ratio, 'nearest');
        end;
        texture = optTexture;
    end
    try
        switch type
            case 'P0'
                p0 = imresize(texture, round(N/n), 'nearest');
                image = p0;                                
            case 'P1'
                width = n;
                height = width;
                p1 = texture(1:height, 1:width);
                image = catTiles(p1, N, type);                
            case 'P2'
                height = round(n/2);
                width = 2*height;
                start_tile = texture(1:height, 1:width);
                tileR180 = imrotate(start_tile, 180);
                p2 = [start_tile; tileR180];
                image = catTiles(p2, N, type);                
            case 'PM'
                height = round(n/2);
                width = 2*height;
                start_tile = texture(1:height, 1:width);
                mirror = flipud(start_tile);
                pm = [start_tile; mirror];
                image = catTiles(pm, N, type);                
            case 'PG'
                height = round(0.5*n);
                width = 2*height;
                start_tile = texture(1:height, 1:width);
                tile = rot90(start_tile, 3);
                glide = flipud(tile);
                pg = [tile, glide];
                image = catTiles(pg, N, type);                
            case 'CM'
                height = round(n/2);
                width = height;
                start_tile = texture(1:height, 1:width);                
                mirror = fliplr(start_tile);
                tile1 = [start_tile, mirror];
                tile2 = [mirror, start_tile];
                cm = [tile1; tile2];
                image = catTiles(cm, N, type);                
            case 'PMM'
                height = round(n/2);
                width = height;
                start_tile = texture(1:height, 1:width);                      
                mirror = fliplr(start_tile);
                pmm = [start_tile, mirror; flipud(start_tile), flipud(mirror)];
                image = catTiles(pmm, N, type);                
            case 'PMG'
                height = round(n/2);
                width = height;
                start_tile = texture(1:height, 1:width);
                pmg = [start_tile, rot90(start_tile, 2); flipud(start_tile), fliplr(start_tile)];
                image = catTiles(pmg, N, type);                
            case 'PGG'
                height = round(n/2);
                width = height;
                start_tile = texture(1:height, 1:width);
                pgg = [start_tile, flipud(start_tile); fliplr(start_tile), rot90(start_tile, 2)];
                image = catTiles(pgg, N, type);                
            case 'CMM'
                height = round(n/4);
                width = 2*height;
                start_tile = texture(1:height, 1:width);                
                tile1 = [start_tile; rot90(start_tile,2)];
                tile2 = flipud(tile1);
                cmm = [tile1, tile2; tile2, tile1];
                image = catTiles(cmm, N, type);                
            case 'P4'
                height = round(n/2);
                width = height;
                start_tile = texture(1:height, 1:width);                
                p4 = [start_tile, rot90(start_tile, 3); rot90(start_tile, 1), rot90(start_tile, 2),];
                image = catTiles(p4, N, type);                
            case 'P4M'
                height = round(n/2);
                width = height;
                start_tile = texture(1:height, 1:width);
                x1 = [0 width 0 0];
                y1 = [0 height height 0];
                mask = poly2mask(x1, y1, width, height);
                tile1 = mask.*start_tile;
                tile2 = rot90(fliplr(tile1), 1);
                tile = max(tile1, tile2);
                p4m = [tile, rot90(tile, 3); rot90(tile, 1), rot90(tile, 2)];
                image = catTiles(p4m, N, type);                
            case 'P4G'
                height = round(n/2);
                width = height;
                start_tile = texture(1:height, 1:width);
                y = [0 width width 0];
                x = [0 0 height 0];
                mask = poly2mask(x, y, height, width);
                tile1 = mask.*start_tile;
                tile2 = rot90(fliplr(tile1), 1);
                tile = max(tile1, tile2);
                p4g = [rot90(tile, 3), rot90(tile, 2); tile, rot90(tile, 1)];
                image = catTiles(p4g, N, type);                
            case 'P3'
                width = round(n/sqrt(3*sqrt(3)));
                height = width;
                start_tile = texture(1:height, 1:width);
                p3 = new_p3(start_tile);
                image = catTiles(p3, N, type);                
            case 'P3M1'
                height = round(n/sqrt(3*sqrt(3)));
                width = height;
                start_tile = texture(1:height, 1:width);                
                p3m1 = new_p3m1(start_tile);                
                image = catTiles(p3m1, N, type);                
            case 'P31M'
                height = round(n/sqrt(sqrt(3)));
                width = height;
                start_tile = texture(1:height, 1:width);                
                %start_tile = ones(height, width);
                p31m = new_p31m(start_tile);
                %ugly trick
                p31m_1 = fliplr(p31m');
                image = catTiles(p31m_1, N, type);                
            case 'P6'
                height = round(n/sqrt(sqrt(3)));
                width = height;
                start_tile = texture(1:height, 1:width);                
                %start_tile = ones(height, width);
                p6 = new_p6(start_tile);
                image = catTiles(p6, N, type);                
            case 'P6M'
                height = round(n/sqrt(sqrt(3)));
                width = height;
                start_tile = texture(1:height, 1:width);
                %start_tile = ones(height, width);
                p6m = new_p6m(start_tile);
                image = catTiles(p6m, N, type);                
            otherwise
                warning('Unexpected Wallpaper Group type. Returning random noise.');
                image = repmat(texture, [ceil(N/n),  ceil(N/n)]);
        end
    catch err
        disp(strcat('new_SymmetricNoise:Error making ', type));
        disp(err.message);
        disp(err.stack(1));
        disp(err.stack(2));
    end
end
function img = catTiles(tile, N, type)
    %disp tile square
    sq = size(tile, 1)*size(tile, 2);
    disp(strcat(type,'_area of tile = ', num2str(sq)));                
    
    %write tile
    imwrite(tile, strcat('~/Documents/MATLAB/tiles/', type, '_tile.jpeg'), 'jpeg');
    dN = 1 + floor(N./size(tile));
    img = repmat(tile, [dN(1), dN(2)]);                
end
%%array of coefficients(DO NOT CHANGE):
%%tile sizes by groups:     tile_in     tile_out         square ratio       width ratio
%P1:                        (n, n)      (n, n)              1                  1
%P2:                        (n, n)      (n, 2n)             2                  0.5
%PM, PG:                    (n, n)      (2n, n)             2                  1
%PMG, PMM, P4, P4M:         (n, n)      (2n, 2n)            4                  0.5
%PGG:                       (n, 2n)     (2n, 2n)            4                  0.5
%CM:                        (n, n)      (n, 2n)             2                  0.25
%P4G:                       (n, 2n)     (4n, 4n)            16                 0.25
%CMM:                       (n, n)      (4n, 4n)            16                 0.25
%P3:                        (n, n)      (3n, n sqrt(3))     3sqrt(3)           1/sqrt(3)
%P31M: s = round(n*2.632):  (s, s)      (3s, s sqrt(3))     3s^2 sqrt(3)/n^2   1/(2.632*sqrt(3))
%P3M1, P6, P6M:             (s, s)      (s, s sqrt(3))      s^2 sqrt(3)/n^2    1/(2.632*sqrt(3))  
   