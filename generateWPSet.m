function generateWPSet

%This script will generate a given number(ngroup) of exemplars of
%given wallpaper groups(defined in groups). Each image will be square(NxN)
%with repeating region of n pixels. Images will be low-pass filtered with

    %% group definitions
    %Groups = {'P1', 'P2', 'PM' ,'PG', 'CM', 'PMM', 'PMG', 'PGG', 'CMM', 'P4', 'P4M', 'P4G', 'P3', 'P3M1', 'P31M', 'P6', 'P6M'};
    %number of images per group
    inGroup = 20;
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
    
    scrambledSet = scrambleGroupSet(rawFreqSet, nScramble);
    
    [psSet, psFreqSet] = psScramble(rawSet);
    
    psScrambledAvgSet = avgPSScrambled(psFreqSet, rawFreqSet);
     
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
    function rawScrambled = scrambleGroupSet(rawFreq, numScramble)
        nGroups = numel(rawFreq);
        rawScrambled = cell(nGroups, 1);
        
        %use average ampliture here
        for i = 1: nGroups           
            rawScrambled{i} = scrambleGroup(rawFreq{i}, numScramble);
        end
    end
    
    %% Scramble group
    function rawScambledGroup = scrambleGroup(inGroup, numScramble)
        rawScambledGroup = cell(numScramble, 1);
        for s = 1:numScramble
            rawScambledGroup{s} = scramble(inGroup);
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
    
    %% replace phase spectrum with random values for given group 
    function outImage = scramble(freqGroup, avgMag)
        if(nargin < 2)
            avgMag = meanMag(freqGroup);
        end;
        randImg = round(255*rand(size(avgMag)));
        randPhase = angle(fft2(double(randImg)));
        
        cmplxIm = avgMag.*exp(1i.*randPhase);
        outImage = ifft2(cmplxIm, 'symmetric');
    end
    
    %% replace each image's amplitude with group average
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
        end
        out = mag/nImages;
    end
    
    %% TODO: analyse group spectrum
    function analyseGroup(rawFreqGroup, scrambledGroup, avgAmp, groupID)
        
        if (nargin < 3)
            avgMag = meanMag(freqGroup);
        end
        
        %plot each image's spectrum, average magnitude and scrambled
        for i = 1:numel(rawFreqGroup)
            freqAnalyser(rawFreqGroup{i}, scrambledGroup{i}, avgAmp, strcat(groupID, 'freq_analysis_', num2str(i)));
        end;
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
