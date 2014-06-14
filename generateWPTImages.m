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
    Groups = {'CMM','P4M'};
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
    timeStr = datestr(now,30);
    timeStr(strfind(timeStr,'T'))='_';
    sPath = strcat(saveStr, timeStr, '/');
    saveFmt = 'png'; %Save fmt/numeration     
    
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
    
    for i = 1:length(Groups)    
        disp(strcat('generating', ' ', Groups{i}));

        group = Groups{i};
        n = round(sqrt(tileArea));
        
        %% generating wallpapers, saving freq. representations
        raw = cellfun(@new_SymmetricNoise,...
            repmat({group},inGroup,1), ...
            repmat({wpSize},inGroup,1),...
            repmat({n},inGroup,1), ...
            'uni',false);
        raw = cellfun(@double,raw,'uni',false);
        rawFreq = cellfun(@fft2,raw,'uni',false);
        
        %% image processing steps
        avgMag = meanMag(rawFreq); % get average magnitude
        avgRaw = cellfun(@spectra,repmat({avgMag},inGroup,1),rawFreq,'uni',false); % replace each image's magnitude with the average 
        filtered = cellfun(@filterImg,avgRaw,repmat({wpSize},inGroup,1),'uni',false); % low-pass filtering + histeq 
        masked = cellfun(@maskImg,filtered,repmat({wpSize},inGroup,1),'uni',false);     % masking the image (final step)
                
        %% making scrambled images
        scrambled_raw = cellfun(@spectra,repmat({avgMag},nScramble,1),'uni',false); % only give spectra only arg, to make randoms
        scrambled_filtered = cellfun(@filterImg,scrambled_raw, repmat({wpSize},inGroup,1),'uni',false);
        scrambled_masked = cellfun(@maskImg,scrambled_filtered,repmat({wpSize},inGroup,1),'uni',false);
        
        %% saving averaged and scrambled images
        groupNumber = mapGroup(group);
        saveStr = arrayfun(@(x) strcat(sPath, 'analysis/plot_',group, '_', num2str(x)),1:inGroup,'uni',false)';
        tempDiff = cellfun(@freqAnalyser, ...
                    repmat({avgMag},inGroup,1),... 
                    avgRaw, filtered,masked, ...
                    scrambled_raw, scrambled_filtered, scrambled_masked, saveStr,'uni',false);
        all_in_one = cellfun(@(x,y,z) cat(2,x(1:wpSize,1:wpSize),y(1:wpSize,1:wpSize),z(1:wpSize,1:wpSize)),...
            raw,avgRaw,filtered,'uni',false);
        for img = 1:inGroup
            if(printAnalysis)
                imwrite(all_in_one{img},  strcat(sPath, 'analysis/steps_',group, '_', num2str(img), '.jpeg'), 'jpeg');
            end;
            patternPath = strcat(sPath,num2str(1000*groupNumber + img), '.', saveFmt);
            saveImg(masked{img},patternPath,saveFmt);
        end
        
        for scr = 1:nScramble
            scramblePath = strcat(sPath,num2str(1000*(groupNumber + 17) + scr), '.', saveFmt);
            saveImg(scrambled_masked{scr},scramblePath,saveFmt);
        end
        if(saveRaw)
            for img = 1:inGroup
                rawPath = strcat(sRawPath,group, '_', num2str(img), '.', saveFmt);
                saveImg(raw{img},rawPath,saveFmt);
            end
        end
        symAveraged(:,i)=[avgRaw;scrambled_raw];
        symFiltered(:,i)= [filtered;scrambled_filtered];
        symMasked(:,i)= [masked;scrambled_masked];
        diffMeans(:,:,i)=cell2mat(tempDiff);
    end
    save([sPath,'analysis/',timeStr,'.mat'],'symAveraged','symFiltered','symMasked','diffMeans','Groups');
end
    
    function saveImg(img,savePath,saveFmt)
        img = uint8(round(img.*255));
        imwrite(img, savePath, saveFmt);
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
    
    %% Portilla - Simoncelli scrambling tool
    % imSet{group}{image} -- cell matrix of raw (unfiltered, uncut) wallpapers
    % out{group}{image} -- scrambled wallpaper in frequency domain  
    function [scrambled, freqScrambled] = psScramble(imSet)

        Nsc = 4; % Number of scales
        Nor = 4; % Number of orientations
        Na = 11;  % Spatial neighborhood is Na x Na coefficients It must be an odd number!

        Niter = 15;	% Number of iterations of synthesis loop
    

        x = numel(imSet);
        scrambled = cell(x, 1);
        freqScrambled = cell(x, 1);
        for i = 1:x
            
            y = numel(imSet{i});
            % two cell arrays to store freq-scrambled and raw-scrambled
            scrambled{i} = cell(y, 1);
            freqScrambled{i} = cell(y, 1);
            
            for j = 1:y
                
                image = imSet{i}{j};
                Nx = size(image, 1);
                Ny = size(image, 2);
                
                %for scrambling algorithm, image dimentions have to be  multiple of 2^(Nsc+2)
                req_mult = 2^(Nsc + 2);
                
                req_Nx = req_mult*floor(Nx/req_mult);
                req_Ny = req_mult*floor(Ny/req_mult);
                % if current wallpaper's size is less than required 
                
                if (req_Nx > Nx || req_Ny > Ny)
                   added_X = cat(2, image, image(:, 1:req_mult));
                   added_XY = cat(1, added_X, added_X(1:req_mult, :));
                   image = added_XY;
                end
                
                Nsx = req_Nx;	
                Nsy = req_Ny;
                mask = zeros(Nsx, Nsy);
                imKeep = zeros(Nsx*Nsy, 1);
                im0 = image(1:Nsx, 1:Nsy);
            
                params = textureAnalysis(im0, Nsc, Nor, Na);

                imKeep(:, 1) = reshape(mask, [Nsx*Nsy, 1]);
                imKeep(:, 2) = reshape(im0, [Nsx*Nsy, 1]);

                res = textureSynthesis(params, [Nsx Nsy], Niter, [], imKeep);
                scrambled{i}{j} = res;
                freqScrambled{i}{j} = fft2(res);
            end
        end
    end
    
    %% replace the PS-scrambled magnitude with group average
    function psScrambledAvgSet = avgPSScrambled(psFreqSet, rawFreqSet)
    
        psScrambledAvgSet = cell(numel(psFreqSet), 1);
        for n = numel(psFreqSet)
            avgMagGroup = meanMag(rawFreqSet{n});
        
            % scrambled wallpapers will have different size and
            % we'll have to crop/cat average magnitude
            % since all images in psScrambledSet{n} have equal size,
            % use the first one to determine it.
        
            Nx = size(psFreqSet{n}{1}, 1);
            Ny = size(psFreqSet{n}{1}, 2);
            diffNx = size(avgMagGroup, 1) - Nx;
            diffNy = size(avgMagGroup, 2) - Ny;
        
            if (diffNx < 0 || diffNy < 0)
                if (diffNx < 0)
                    added_X = cat(2, avgMagGroup, avgMagGroup(:, 1:abs(diffNx)));
                    added_XY = cat(1, added_X, added_X(1:abs(diffNy), :));
                    avgMagGroup = added_XY;
                end
            end
            avgMagGroup = avgMagGroup(1:Nx, 1:Ny);
            psScrambledAvgSet{n} = meanGroup(psFreqSet{n}, avgMagGroup);
        end
    end


