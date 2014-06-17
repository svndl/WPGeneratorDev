function outputSet = replaceMagnitude(inputSet, newMag, domain)

% For every image in inputSet replaces the magnitude with newMag
% if newMag is not specified, average magnitude will be used for repacing
% inputSet and newMag used as if they are defined in same domain
% outSet is  


    try
        if (nargin < 2)
            newMag = averageMagnitude(inputSet, domain);
        end;
        
        nImages = length(inputSet);
        outputSet = cell(nImages, 1);
        
        switch domain
            case 'FREQ'
                for n = 1:nImages
                    outputSet{n} = ifft2(newMag.*exp(1i.*angle(inputSet{n})), 'symmetric');
                end
            case 'SPACE'
                for n = 1:nImages
                    outputSet{n} = ifft2(newMag.*exp(1i.*angle(fft2(inputSet{n}))), 'symmetric');
                end
            otherwise
                disp('WARNING: domain is not specified, assuming FREQ');
                outputSet = replaceMagnitude(inputSet, newMag, 'FREQ');
        end
    catch err
        disp('replaceMagnitude:Error  ');
        disp(err.message);
        disp(err.stack(1));
        disp(err.stack(2));
        outputSet = 0;
    end
end
