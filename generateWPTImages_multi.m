function generateWPTImages_multi(groupNames, n_exemplars, wp_sizes, mask_image, save_fmt, run_par, out_dir)
    if nargin < 7
        out_dir = '~/Desktop/WPset';
    else
    end
    if nargin < 6
        run_par = true;
    else
    end
    if nargin < 5
        save_fmt = 'PGM'; %Save fmt/numeration     
    else
    end
    if nargin < 4
        mask_image = true;
    else
    end
    if nargin < 3
        % image size and tile size
        wp_sizes = [600,100];
    else
    end
    if nargin < 2
        n_exemplars = 10;
    else
    end   
    if nargin < 1
        groupNames = {'P1', 'P2', 'PM' ,'PG', 'CM', 'PMM', 'PMG', 'PGG', 'CMM', 'P4', 'P4M', 'P4G', 'P3', 'P3M1', 'P31M', 'P6', 'P6M'};
    else
    end
    
    % if not cell, wrap
    if ~iscell(groupNames)
        groupNames = {groupNames};
    else
    end
    
    %% check if matlab pool is open, close, then open
    if run_par
        if ~isempty(gcp('nocreate'))
            delete(gcp('nocreate'));
        else
        end
        parpool(4)
        n_workers = 4;
    else
        n_workers = 0;
    end
    
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
    %number of images per group
    inGroup = n_exemplars;
    
    %% image parameters
    if numel(wp_sizes) < 2
        % assume a ratio of 6
        wp_sizes(2) = wp_sizes(1)/6;
    else
    end
    %area of tile that will be preserved across groups
    tileArea = wp_sizes(2)*wp_sizes(2);    
    
    %% define number of scrambled images per group
    nScramble = n_exemplars;    

    %% Average magnitude within the each group
    %%save parameters
    saveStr = sprintf('%s/dev',out_dir);
    timeStr = datestr(now,30);
    timeStr(strfind(timeStr,'T'))='_';
    sPath = strcat(saveStr, timeStr, '/');
    
    
    
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
    printAnalysis = false;
    
    try
        mkdir(sPath);
        if(saveRaw)
            mkdir(sRawPath);
        end
        if(printAnalysis)
            mkdir(sAnalysisPath)
        end
    catch err
        error('MATLAB:generateWPSet:mkdir', sPath);
    end
    
    %% Generating WPs and scrambling
    
    parfor (i = 1:length(groupNames), n_workers)    
        disp(strcat('generating', ' ', groupNames{i}));

        group = groupNames{i};
        n = round(sqrt(tileArea));
        
        %% generating wallpapers, saving freq. representations
        raw = cellfun(@new_SymmetricNoise,...
            repmat({group},inGroup,1), ...
            repmat({wp_sizes(1)},inGroup,1),...
            repmat({n},inGroup,1), ...
            'uni',false);
        raw = cellfun(@double,raw,'uni',false);
        rawFreq = cellfun(@fft2,raw,'uni',false);
        
        %% image processing steps
        avgMag = meanMag(rawFreq); % get average magnitude
        avgRaw = cellfun(@spectra,repmat({avgMag},inGroup,1),rawFreq,'uni',false); % replace each image's magnitude with the average 
        filtered = cellfun(@filterImg,avgRaw,repmat({wp_sizes(1)},inGroup,1),'uni',false); % low-pass filtering + histeq 
        masked = cellfun(@maskImg,filtered,repmat({wp_sizes(1)},inGroup,1),'uni',false);     % masking the image (final step)
                
        %% making scrambled images
        scrambled_raw = cellfun(@spectra,repmat({avgMag},nScramble,1),'uni',false); % only give spectra only arg, to make randoms
        scrambled_filtered = cellfun(@filterImg,scrambled_raw, repmat({wp_sizes(1)},inGroup,1),'uni',false);
        scrambled_masked = cellfun(@maskImg,scrambled_filtered,repmat({wp_sizes(1)},inGroup,1),'uni',false);
        
        %% saving averaged and scrambled images
        groupNumber = mapGroup(group);
        if printAnalysis
            saveStr = arrayfun(@(x) strcat(sPath, 'analysis/plot_',group, '_', num2str(x)),1:inGroup,'uni',false)';
        else
            saveStr = arrayfun(@(x) false,1:inGroup,'uni',false)';
        end
        tempDiff = cellfun(@freqAnalyser, ...
                    repmat({avgMag},inGroup,1),... 
                    avgRaw, filtered,masked, ...
                    scrambled_raw, scrambled_filtered, scrambled_masked, saveStr,'uni',false);
        all_in_one = cellfun(@(x,y,z) cat(2,x(1:wp_sizes(1),1:wp_sizes(1)),y(1:wp_sizes(1),1:wp_sizes(1)),z(1:wp_sizes(1),1:wp_sizes(1))),...
            raw,avgRaw,filtered,'uni',false);
        for img = 1:inGroup
            if(printAnalysis)
                imwrite(all_in_one{img},  strcat(sPath, 'analysis/steps_',group, '_', num2str(img), '.jpeg'), 'jpeg');
            end;
            patternPath = strcat(sPath,num2str(1000*groupNumber + img), '.', save_fmt);
            if mask_image
                saveImg(masked{img},patternPath,save_fmt);
            else
                saveImg(filtered{img},patternPath,save_fmt);
            end
               
        end
        
        for scr = 1:nScramble
            scramblePath = strcat(sPath,num2str(1000*(groupNumber + 17) + scr), '.', save_fmt);
            if mask_image
                saveImg(scrambled_masked{scr},scramblePath, save_fmt);
            else
                saveImg(scrambled_filtered{img},scramblePath, save_fmt);
            end
        end
        if(saveRaw)
            for img = 1:inGroup
                rawPath = strcat(sRawPath,group, '_', num2str(img), '.', save_fmt);
                saveImg(raw{img},rawPath,save_fmt);
            end
        end
        symName = groupNames{i};
        symAveraged = [avgRaw;scrambled_raw];
        symFiltered = [filtered;scrambled_filtered];
        symMasked = [masked;scrambled_masked];
        diffMeans = cell2mat(tempDiff);
        saveMat(symName,symAveraged,symFiltered,symMasked,diffMeans,sPath,timeStr);
    end
    delete(gcp('nocreate'));
end

    function saveMat(name,averaged,filtered,masked,diffmeans,sPath,timeStr)
        group.Name = name;
        group.Averaged = averaged;
        group.Filtered = filtered;
        group.Masked = masked;
        group.diffMeans = diffmeans;
        save([sPath,'/',group.Name,'_',timeStr,'.mat'],'group','-v7.3');
    end 

    function saveImg(img,savePath,save_fmt)
        img = uint8(round(img.*255));
        imwrite(img, savePath, save_fmt);
    end 

    %% Filter/mask every image
    function outImg = filterImg(inImg, N)        
        % Make filter intensity adaptive (600 is empirical number)
        sigma = N/600;
        lowpass = fspecial('gaussian', [9 9], sigma);
    
        % filter
        image = imfilter(inImg, lowpass);
        
        % histeq
        image = histeq(image);
                
        % normalize
        image = (image)./range(image(:)); %scale to unit range
        image = image - mean(image(:)); %bring mean luminance to zero		
        image = image/max(abs(image(:))); %Scale so max signed value is 1
        image = 125*image+127; % Scale into 2-252 range
        image = image./255;
        
        outImg = image;
    end
    
    %% apply mask
    function outImg = maskImg(inImg, N)
        %define mask(circle)
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
    
    %% save group
    function writeGroup(path, type, data, save_fmt, extra)
        if(nargin < 5)
            extra = '';
        end;
        nImages = length(data);
        for n = 1:nImages
            filename = strcat(type, num2str(n), extra, '.', save_fmt);              
            imwrite(data{n}, strcat(path, filename), save_fmt);
        end
    end
    
    %% replace spectra
    function outImage = spectra(avgMag,imFreq)
        if(nargin < 2) % if no image frequency input, make random image and get the frequency
            randImg = rand(size(avgMag));
            imFreq = fft2(double(randImg));
        end
        cmplxIm = avgMag.*exp(1i.*angle(imFreq));
        outImage = ifft2(cmplxIm, 'symmetric');
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
