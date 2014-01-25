function generateWPSet

%This script will generate a given number(ngroup) of exemplars of
%given wallpaper groups(defined in groups). Each image will be square(NxN)
%with repeating region of n pixels. Images will be low-pass filtered with
%lowpass

    %%group definitions
    %groups = {'P1', 'P2', 'PM' ,'PG', 'CM', 'PMM', 'PMG','PGG', 'CMM', 'P4', 'P4M', 'P4G','P3', 'P3M1', 'P31M', 'P6', 'P6M'};
    %groups = {'P3M1', 'P31M', 'P6', 'P6M'};
    %groups = {'P1', 'P2', 'PM' ,'PG', 'CM', 'PMM', 'PMG','PGG', 'CMM', 'P4', 'P4M', 'P4G','P3', 'P6', 'P6M'};
    %groups = {'P6', 'P6M'};
    groups = {'P3', 'P1', 'P4', 'P2'};
    %number of images per group
    %groups = {'P3'};
    ingroup = 5;
    ngroups = length(groups);
    
    %%image parameters
    %image size
    N = 1024;
    
    %area of tile that will be preserved across groups
    n = 160;    
    
    %%Define your texture here
    %myTexture = ...
    
    %%save parameters
    saveStr = '~/Documents/WPSet/dev/';
    sPath = strcat(saveStr, datestr(clock), '/');
    
    %% define number of scrambled images per group
    nScramble = 5;
    
    %%store the raw images
    sRaw = true;
    sRawPath = strcat(sPath, 'raw/');
    
    %%cell array to store raw images per group
    raw = cell(ingroup, 1);
    
    %cell array to store ffts of images per group
    rawFreq = cell(ingroup, 1);
    
    %cell array to store scrambled
    rawScambled = cell(nScramble, 1);
    try
        mkdir(sPath);
        if(sRaw)
            mkdir(sRawPath);
        end;
    catch err
        error('MATLAB:generateWPSet:mkdir', sPath);
    end;
    
    %Generating WPs and scrambling
    for i = 1:ngroups
        
        group = groups{i};
        disp(strcat('generating ', group));
        
        %%generating wallpapers, saving freq. representations
        for j = 1:ingroup
            raw{j} = new_SymmetricNoise(group, N, n);
            rawFreq{j} = fft2(double(raw{j})); 
        end;
        
        %save average magnitude
        avgMag = meanMag(rawFreq);
        imwrite(avgMag, strcat(sPath, group, '_magnitude.jpeg'), 'jpeg');

        %averaging images (replace each image's magnitude with the average) 
        averaged = filterGroup(meanGroup(rawFreq, avgMag), N);
        filtered = filterGroup(raw, N);
        
        %saving images
        writeGroup(sPath, group, filtered);
        writeGroup(sPath, group, averaged, '_averaged');

        %%scrambling 
        for s = 1:nScramble
            rawScambled{s} = scramble(rawFreq, avgMag);
        end
        scrambled = filterGroup(rawScambled, N);
        writeGroup(sPath, group, scrambled, '_scrambled');
        
        if (sRaw)
           writeGroup(sRawPath, group, raw);
           writeGroup(sRawPath, group, rawScambled, 'scrambled');
        end
    end
end
    
    %%apply mask/filter/histeq to group 
    function out = filterGroup(inGroup, N)
        num = length(inGroup);
        out = cell(num, 1);
        for n = 1:num
            out{n} = filterImg(inGroup{n}, N);
        end;
    end
    
    %%Filter/mask every image
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
    
    %%save group
    function writeGroup(path, type, data, extra)
        if(nargin < 4)
            extra = '';
        end;
        nImages = length(data);
        for n = 1:nImages
            filename = strcat(type, '_', num2str(n), extra, '.jpeg');              
            imwrite(data{n}, strcat(path, filename), 'jpeg');
        end
    end
    
    %%replace spectra
    
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
    
    %returns average mag of the group
    function out = meanMag(freqGroup)
        nImages = length(freqGroup);
        mag = 0;
        for n = 1:nImages
            mag = mag + abs(freqGroup{n});
        end;
        out = mag/nImages;
    end
    
    function analyseGroup(rawFreqGroup, scrambledGroup, avgAmp)
        if (nargin < 3)
            avgMag = meanMag(freqGroup);
        end
        
        %plot each image's spectrum, average magnitude and scrambled
    end