function generateWPTImages_multi
    %% check if matlab pool is open, close, then open
    if matlabpool('size') > 0 
        matlabpool close
    else
    end
    matlabpool open 4
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
    Groups = keySet;
    %number of images per group
    inGroup = 20;
    
    %% image parameters
    %image size
    wpSize = 600;
    %area of tile that will be preserved across groups
    tileArea = 100*100;    
    
    %% define number of scrambled images per group
    nScramble = 20;    

    %% Average magnitude within the each group
    %%save parameters
    saveStr = '~/Documents/WPSet/dev/';
    sPath = strcat(saveStr, datestr(clock), '/');
    saveFmt = 'pgm'; %Save fmt/numeration     
    
    %% Handling raw images 
    saveRaw = false;
    sRawPath = strcat(sPath, 'raw/');
    sAnalysisPath = strcat(sPath, 'analysis/');
    
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
        if(printAnalysis)
            mkdir(sAnalysisPath)
        end
    catch err
        error('MATLAB:generateWPSet:mkdir', sPath);
    end;
    
    %% Generating WPs and scrambling
    
    parfor i = 1:length(Groups)    
        disp(strcat('generating', ' ', Groups{i}));
        
        group = Groups{i};
        n = round(sqrt(tileArea));
        %%generating wallpapers, saving freq. representations
        raw = cellfun(@new_SymmetricNoise,...
            repmat({group},inGroup,1), ...
            repmat({wpSize},inGroup,1),...
            repmat({n},inGroup,1), ...
            'uni',false);
        raw = cellfun(@double,raw,'uni',false);
        rawFreq = cellfun(@fft2,raw,'uni',false); 
        %averaging images (replace each image's magnitude with the average) 
        avgMag = meanMag(rawFreq);
        avgRaw = meanGroup(rawFreq, avgMag);
        filtered = filterGroup(avgRaw, wpSize);
        masked = maskGroup(filtered, wpSize);
        
        %% saving averaged and scrambled images
        groupNumber = mapGroup(group);
        for img = 1:inGroup
            if(printAnalysis) 
                saveStr = strcat(sPath, 'analysis/',group, '_', num2str(img), 'analysis');
                freqAnalyser(avgMag, rawFreq{img}, filtered{img}, saveStr);
                all_in_one = cat(2, raw{img}(1:wpSize, 1:wpSize), avgRaw{img}(1:wpSize, 1:wpSize), filtered{img}(1:wpSize, 1:wpSize));
                imwrite(all_in_one,  strcat(sPath, 'analysis/',group, '_', num2str(img), '.jpeg'), 'jpeg');
            end;
            wallpaperName = strcat(num2str(1000*groupNumber + img), '.', saveFmt);
            imwrite(filtered{img}, strcat(sPath, wallpaperName), saveFmt);
        end
        
        %% making scrambled images
        rawScrambled = cellfun(@scramble,repmat({rawFreq},nScramble,1),repmat({avgMag},nScramble,1),'uni',false);
        scrambled = filterGroup(rawScrambled, wpSize);

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
    matlabpool close
end


    %% Filter/histeq group 
    function out = filterGroup(inGroup, N)
        num = length(inGroup);
        out = cell(num, 1);
        for n = 1:num
            out{n} = filterImg(inGroup{n}, N);
        end;
    end
    
    %% Mask group
    function out = maskGroup(inGroup, N)
        num = length(inGroup);
        out = cell(num, 1);
        for n = 1:num
            out{n} = maskImg(inGroup{n}, N);
        end
    end
    %% Filter/histeq image
    function outImg = filterImg(inImg, N)        
        % Make filter intensity adaptive (600 is empirical number)
        sigma = N/600;
        lowpass = fspecial('gaussian', [9 9], sigma);
    
        %%filter
        image = imfilter(inImg, lowpass);
        
        %normalize
        image = image - min(image(:));
        image = image./max(image(:));
        %%histeq
        image_hn = histeq(image);
        outImg = image_hn;
        %prepare for PowerDiva (remove 1s and 0s), keep centered around 0.5        
        outImg(outImg<1/255) = 1/255; 
    end
    
    %% Mask image
    function outImg = maskImg(inImg, N)
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
        outImg = inImg(1:size(mask, 1), 1:size(mask, 2));
        outImg(mask==0)=.5;
    end
    
    %% Save group
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
