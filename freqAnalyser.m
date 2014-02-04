function freqAnalyser(meanAmp, currAmp, histeqAmp, saveString)
    h = figure('visible','off');
    
    oriSpecAxH = subplot(2, 1, 1);
    sfSpecAxH = subplot(2, 1, 2);
    
    
    mOriPow = meanOriPow(meanAmp);
    cOriPow = meanOriPow(currAmp);
    cHisteqPow = meanOriPow(fft2(histeqAmp));

    yScale = linspace(-180, 180, 361);
    %% disp Orientation energy
    semilogy(oriSpecAxH, yScale, [mOriPow, cOriPow, cHisteqPow]);
    title(oriSpecAxH, 'Orientation Energy');
    legend(oriSpecAxH, 'mean', 'current', 'mean histeq');
    xlabel(oriSpecAxH, 'Orientation (deg)');
    
    %% disp spatial energy
    mSpatialPow = meanSpatialPow(meanAmp);
    cSpatialPow = meanSpatialPow(currAmp);
    hSpatialPower = meanSpatialPow(fft2(histeqAmp));
    
    frqs = linspace(0,size(meanAmp, 1)/2,101);
    loglog(sfSpecAxH, frqs, [mSpatialPow, cSpatialPow, hSpatialPower]);
    title(sfSpecAxH, 'Spatial Frequency Power');
    legend(sfSpecAxH, 'mean', 'current', 'mean histeq');        
    xlabel(sfSpecAxH, 'Cycles per Image');
    
    saveas(h, saveString, 'jpeg');
    close(h);
end
    
function out = meanOriPow(inImage)
    [nx, ny] = size(inImage);
    [x, y] = ndgrid(linspace(-1, 1, nx),linspace(-1, 1, ny));
    [th, r] = cart2pol(x, y);
    
    th = round( (180/pi)*th);
    th = th(:);
    out = zeros(361, 1);
    idx = 1;
        
    f = fftshift(inImage);
    th(r < .1) = 255;
    th(r > 1) = 255; 
    for i = -180:180,            
        %normFac = max(r(th==i));
        out(idx) = mean(abs(f(th == i)));
        idx = idx + 1;
    end
end
function out = meanSpatialPow(inImage)
    [nx, ny] = size(inImage);
    [x, y] = ndgrid(linspace(-1, 1, nx), linspace(-1, 1, ny));
    [th, r] = cart2pol(x,y);
    
    r = round(100*r);
    out = zeros(101, 1);
    idx = 1;
    f = fftshift(inImage);
    for i = 0:100,
        out(idx) = mean(abs(f(r == i)));
        idx = idx + 1;
    end
end
