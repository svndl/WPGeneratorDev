function outTile =  filterTile(inTile, filterIntensity)
    %outTile = generate_ftile(size(inTile, 1), size(inTile, 2));

    mu = 0.5;
    nx = size(inTile, 1);
    ny = size(inTile, 2);
    %% make adaptive filtering
    
    sigma_x = 10*filterIntensity/nx;
    sigma_y = 10*filterIntensity/ny;
    x = linspace(0, 1, nx);
    y = linspace(0, 1, ny);
    
    gx = exp(-((x - mu).^2)/(2*sigma_x^2))/(sigma_x*sqrt(2*pi));
    gy = exp(-((y - mu).^2)/(2*sigma_y^2))/(sigma_y*sqrt(2*pi));
    
    
    gauss2 = gx'*gy;
    gauss2 = gauss2 - min(min(gauss2));
    gauss2 = gauss2/max(max(gauss2));
    gauss2 = gauss2*5;
    filtered = abs(ifft2(fft2(inTile).*gauss2));
    
    %normalize tile

    outTile = filtered - min(min(filtered));
    outTile = outTile/max(max(outTile));
    %outTile = histeq(outTile);
    
end
