function generateWPTImages
    %% define group to index mapping
    keySet = {'P1', 'P2', 'PM' ,'PG', 'CM', 'PMM', 'PMG', 'PGG', 'CMM', 'P4', 'P4M', 'P4G', 'P3', 'P3M1', 'P31M', 'P6', 'P6M'};
    valueSet = 101:1:117;
    mapGroup = containers.Map(keySet, valueSet);
    
    hexLattice = {'P3', 'P3M1', 'P31M', 'P6', 'P6M'};
    sqrLattice = {'P4', 'P4M', 'P4G'};
    recLattice = {'PM', 'PMM', 'PMG', 'PGG', 'PG'};
    rhoLattice = {'CM', 'CMM'};
    obqLattice = {'P1', 'P2'};
    %% define groups to be generated
    Groups = hexLattice;
    %number of images per group
    inGroup = 5;
    
    %% image parameters
    %image size
    wpSize = 600;
    %area of tile that will be preserved across groups
    tileArea = 100*100;    
    
    %% define number of scrambled images per group
    nScramble = 2;    

    %% Average magnitude within the each group
    %%save parameters
    saveStr = '~/Documents/WPSet/dev/';
    sPath = strcat(saveStr, datestr(clock), '/');
    saveFmt = 'pgm'; %Save fmt/numeration     
    
    
    %% Handling raw images 
    saveRaw = false;
    sRawPath = strcat(sPath, 'raw/');
    
    %cell array to store raw images per group
    raw = cell(inGroup, 1);
    
    %cell array to store ffts of images per group
    rawFreq = cell(inGroup, 1);
    
    %cell array to store scrambled
    rawScambled = cell(nScramble, 1);
    printAnalysis = true;
    
    try
        mkdir(sPath);
        if(saveRaw)
            mkdir(sRawPath);
        end;
    catch err
        error('MATLAB:generateWPSet:mkdir', sPath);
    end;
    
    %% Generating WPs and scrambling
    
    for i = 1:length(Groups)    
        group = Groups{i};
        disp(strcat('generating', ' ', group));
        n = round(sqrt(tileArea));
        
        %%generating wallpapers, saving freq. representations
        for j = 1:inGroup
            raw{j} = new_SymmetricNoise(group, wpSize, n);
            rawFreq{j} = fft2(double(raw{j})); 
        end;
        
        %averaging images (replace each image's magnitude with the average) 
        avgMag = meanMag(rawFreq);
        avgRaw = meanGroup(rawFreq, avgMag);
        averaged = filterGroup(avgRaw, wpSize); 
        
 
        %% saving averaged and scrambled images
        groupNumber = mapGroup(group);
        for img = 1:inGroup
            if(printAnalysis) 
                saveStr = strcat(sPath, group, '_', num2str(img), 'analysis');
                freqAnalyser(avgMag, rawFreq{img}, averaged{img}, saveStr);
                all_in_one = cat(2, raw{img}(1:wpSize, 1:wpSize), avgRaw{img}(1:wpSize, 1:wpSize), averaged{img});
                imwrite(all_in_one,  strcat(sPath, group, '_', num2str(img), '.jpeg'), 'jpeg');
            end;
            wallpaperName = strcat(num2str(1000*groupNumber + img), '.', saveFmt);
            imwrite(averaged{img}, strcat(sPath, wallpaperName), saveFmt);
        end
        
        %% making scrambled images
        for s = 1:nScramble
            rawScambled{s} = scramble(rawFreq, avgMag);
        end
        scrambled = filterGroup(rawScambled, wpSize);

        
        for scr = 1:nScramble
            scrambleName = strcat(num2str(1000*(groupNumber + 17) + scr), '.', saveFmt);
            imwrite(scrambled{scr}, strcat(sPath, scrambleName), saveFmt);
        end
        if(saveRaw)
            for img = 1:inGroup
                wallpaperName = strcat(group, '_', num2str(img), '.', saveFmt);
                imwrite(raw{img}, strcat(sRawPath, wallpaperName), saveFmt);
            end
        end
        allMeaned(:,i)=avgRaw;
        allAveraged(:,i)= averaged;
        allScrambled(:,i)= scrambled;
        allRaw(:,i)=raw;
    end
    save([sPath,datestr(clock),'.mat'],'allAveraged','allScrambled','allRaw','allMeaned','Groups');
end


    %% apply mask/filter/histeq to group 
    function out = filterGroup(inGroup, N)
        num = length(inGroup);
        out = cell(num, 1);
        for n = 1:num
            out{n} = filterImg(inGroup{n}, N);
        end;
    end
    
    %% Filter/mask every image
    function outImg = filterImg(inImg, N)
    %filter parameters(whole image)
        
        % Make filter intensity adaptive (600 is empirical number)
        sigma = N/600;
        lowpass = fspecial('gaussian', [9 9], sigma);
        %%define mask(circle)
        r = 0.5*N;
        X = -0.5*N:0.5*N - 1;   
        X = repmat(X, [N, 1]);
        Y = X';
        D = sqrt(X.^2 + Y.^2);
        D = D./r;
        D(D < 1) = 0;
        D(D > 1) = 1;
        mask = 1 - D;

    
        %%filter
        image = imfilter(inImg, lowpass);
        
        %normalize
        image = image - min(image(:));
        image = image./max(image(:));
        %%histeq
        image_hn = histeq(image);
       
        %%apply mask
        %masked = 1 - image_hn(1:size(mask, 1), 1:size(mask, 2)).*mask;
        masked = image_hn(1:size(mask, 1), 1:size(mask, 2));
        masked(mask==0) = 1;
        
        %prepare for PowerDiva (remove 0s), keep centered around 0.5
        outImg = abs(masked) + 1/255; 
        grayMask = 0.5*(1 - mask);
        outImg = outImg - grayMask;
    end
    
    %% save group
    function writeGroup(path, type, data, saveFmt, extra)
        if(nargin < 5)
            extra = '';
        end;
        nImages = length(data);
        for n = 1:nImages
            filename = strcat(type, num2str(n), extra, '.', saveFmt);              
            imwrite(data{n}, strcat(path, filename), saveFmt);
        end
    end
    
    %% replace spectra
    
    function outImage = scramble(freqGroup, avgMag)
        if(nargin < 2)
            avgMag = meanMag(freqGroup);
        end;
        randImg = rand(size(avgMag));
        randPhase = angle(fft2(double(randImg)));
        
        cmplxIm = avgMag.*exp(1i.*randPhase);
        outImage = ifft2(cmplxIm, 'symmetric');
    end
    
    function out = meanGroup(freqGroup, avgMag)
        if (nargin < 2)
            avgMag = meanMag(freqGroup);
        end;
        nImages = length(freqGroup);
        out = cell(nImages, 1);
        for n = 1:nImages
            out{n} = ifft2(avgMag.*exp(1i.*angle(freqGroup{n})), 'symmetric');
        end
    end
    
    %% returns average mag of the group
    function out = meanMag(freqGroup)
        nImages = length(freqGroup);
        mag = [];
        for n = 1:nImages
            mag(:,:,n) = abs(freqGroup{n});
        end;
        out = median(mag,3);
    end
