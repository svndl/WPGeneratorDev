function generateWPSet

%This script will generate a given number(ngroup) of exemplars of
%given wallpaper groups(defined in groups). Each image will be square(NxN)
%with repeating region of n pixels. Images will be low-pass filtered with

    %% group definitions
    Groups = {'P1', 'P2', 'PM' ,'PG', 'CM', 'PMM', 'PMG', 'PGG', 'CMM', 'P4', 'P4M', 'P4G', 'P3', 'P3M1', 'P31M', 'P6', 'P6M'};
    %number of images per group
    inGroup = 5;
    
    %% image parameters
    %image size
    wpSize = 1024;
    %area of tile that will be preserved across groups
    tileArea = 160*160;    
    
    %% define number of scrambled images per group
    nScramble = 5;
    
    %% remember current rand settings, switch to shuffle
    scurr =  rng;
    rng('shuffle');
    
    %% Generate raw and scrambled set of groups 
    [rawSet, rawFreqSet] = generateGroupSet(Groups, inGroup, wpSize, tileArea);
    
    %scrambledSet = scrambleGroupSet(rawFreqSet, nScramble);
    
    %% average magnitude
    
    %% save averaged ans scrambled sets
    
    %% switch back to prev rand settings
    rng(scurr);
    
    %% Average magnitude within the each group
    %%save parameters
    saveStr = '~/Documents/WPSet/dev/';
    sPath = strcat(saveStr, datestr(clock), '/');
    saveMode = 'jpeg'; %Save fmt/numeration     
    
    
    %% Handling raw images 
    sRaw = true;
    sRawPath = strcat(sPath, 'raw/');
    
    try
        mkdir(sPath);
        if(sRaw)
            mkdir(sRawPath);
        end;
    catch err
        error('MATLAB:generateWPSet:mkdir', sPath);
    end;
    
    %% Saving WPs 
    for i = 1:length(Groups)
        res = filterGroup(rawSet{i}, wpSize);
        writeGroup(sPath, Groups{i}, res, saveMode);
        if (sRaw)
            writeGroup(sRawPath, Groups{i}, rawSet{i}, saveMode);
        end
    end    
end

    %% 
    
    function [rawSet, rawFreqSet] = generateGroupSet(groups, inGroup, imgSize, tileArea)
        %num of groups
        nGroups = length(groups);
        %cell array to store each group's set
        rawSet = cell(nGroups, 1);
        rawFreqSet = cell(nGroups, 1);
        
        for i = 1:nGroups
            [rawSet{i}, rawFreqSet{i}] = generateGroup(groups{i}, inGroup, imgSize, tileArea);
        end
    end
    
    %% generate a given number of exemplars within the typeGroup
    function [rawGroup, rawFreqGroup] = generateGroup(typeGroup, inGroup, imgSize, tileArea)
        
        rawGroup = cell(inGroup, 1);
        rawFreqGroup = cell(inGroup, 1);
        n = round(sqrt(tileArea));
        
        disp(strcat('generating ', typeGroup));
        for j = 1:inGroup
            rawGroup{j} = new_SymmetricNoise(typeGroup, imgSize, n);
            rawFreqGroup{j} = fft2(double(rawGroup{j})); 
        end;
    end
   
    %% Scramble wallpaper set
    function rawScambled = scrambleWallpapers(rawFreq, numScramble)
    end
    
    %% Scramble group
    function rawScambledGroup = scrambleGroup(inGroup, numScramble)
        for s = 1:nScramble
            rawScambledGroup{s} = scramble(rawFreq, avgMag);
        end
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
    
        lowpass = fspecial('gaussian', [9 9], 1.5);
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
        image = image - min(image(:));
        image = image./max(image(:));
        
        %%histeq
        image_hn = histeq(image);
        %%apply mask
        outImg = 1 - image_hn(1:size(mask, 1), 1:size(mask, 2)).*mask;
    end
    
    %% save group
    function writeGroup(path, type, data, saveFmt, extra)
        if(nargin < 5)
            extra = '';
        end;
        nImages = length(data);
        for n = 1:nImages
            filename = strcat(type, '_', num2str(n), extra, '.', saveFmt);              
            imwrite(data{n}, strcat(path, filename), saveFmt);
        end
    end
    
    %% replace spectra
    
    function outImage = scramble(freqGroup, avgMag)
        if(nargin < 2)
            avgMag = meanMag(freqGroup);
        end;
        randImg = round(255*rand(size(avgMag)));
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
        mag = 0;
        for n = 1:nImages
            mag = mag + abs(freqGroup{n});
        end;
        out = mag/nImages;
    end
    
    %% TODO: analyse group spectrum
    function analyseGroup(rawFreqGroup, scrambledGroup, avgAmp)
        if (nargin < 3)
            avgMag = meanMag(freqGroup);
        end
        
        %plot each image's spectrum, average magnitude and scrambled
    end
