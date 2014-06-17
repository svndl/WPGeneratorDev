function generateWPSet

%This script will generate a given number(ngroup) of exemplars of
%given wallpaper groups(defined in groups). Each image will be square(NxN)
%with repeating region of n pixels. Images will be low-pass filtered with

    %% group definitions
    %Groups = {'P1', 'P2', 'PM' ,'PG', 'CM', 'PMM', 'PMG', 'PGG', 'CMM', 'P4', 'P4M', 'P4G', 'P3', 'P3M1', 'P31M', 'P6', 'P6M'};
    %number of images per group
    inGroup = 3;
    Groups = {'P2', 'P3', 'P4', 'P6'};
    
    %% image parameters
    %image size
    wpSize = 800;
    %area of tile that will be preserved across groups
    tileArea = 160*160;    
    
    %% define number of scrambled images per group
    nScramble = 5;
    
    %% Remember current rand generation settings and switch to true random
    
    scurr =  rng;
    rng('shuffle');

    %% remember current rand settings, switch to shuffle
    scurr =  rng;
    rng('shuffle');
    
    %% Generate raw and scrambled set of groups 
    [rawSet, rawFreqSet] = generateGroupSet(Groups, inGroup, wpSize, tileArea);
    
    scrambledSet = scrambleGroup(rawFreqSet, nScramble);
    
    %[psSet, psFreqSet] = psScramble(rawSet);
    
    %psScrambledAvgSet = avgPSScrambled(psFreqSet, rawFreqSet);
     
    %% save averaged ans scrambled sets
    

    %% restore rand settings
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
    sScrambled = true;
    sScrambledPath = strcat(sPath, 'scrambled/');
    
    %% Create directories

    try
        mkdir(sPath);
        if(sRaw)
            mkdir(sRawPath);
        end
        if(sScrambled)
            mkdir(sScrambledPath);
        end
    catch err
        error('MATLAB:generateWPSet:mkdir', sPath);
    end;
    
    %% Saving WPs 
    for i = 1:length(Groups)
        resSet = filterGroup(rawSet{i}, wpSize);
        writeGroup(sPath, Groups{i}, resSet, saveMode);
        if (sRaw)
            writeGroup(sRawPath, Groups{i}, rawSet{i}, saveMode);
            writeGroup(sRawPath, Groups{i}, psSet{i}, saveMode, '_PSscrambled');
            writeGroup(sRawPath, Groups{i}, psScrambledAvgSet{i}, saveMode, '_avgPSscrambled');
        end;
        if(sScrambled)
            resScrambled = filterGroup(scrambledSet{i}, wpSize);
            resPSScrambled = filterGroup(psSet{i}, wpSize);
            resPSAvgScrambled = filterGroup(psScrambledAvgSet{i}, wpSize);
            writeGroup(sScrambledPath, Groups{i}, resScrambled, saveMode, '_scrambled');
            writeGroup(sScrambledPath, Groups{i}, resPSAvgScrambled, saveMode, '_psAvgScrambled');
            writeGroup(sScrambledPath, Groups{i}, resPSScrambled, saveMode, '_psScrambled');            
        end
    end
end

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
    function rawScrambled = scrambleGroup(rawFreq, numScramble)
        nGroups = length(rawFreq);
        rawScrambled = cell(nGroups, 1);
        
        %use average ampliture here
        for i = 1:nGroups 
            avgMag = averageMagnitude(rawFreq{i}, 'FREQ');
            rawScrambled{i} = scramblePhase(avgMag, numScramble);
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
    