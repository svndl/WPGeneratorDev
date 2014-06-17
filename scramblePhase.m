function outputSet = scramblePhase(inputImgMagnitude, nScramble)
% for each inputImgMagnitude generates nScramble images with inputImage magnitude and random phase
% input magnitude is expected in frequency domain, output set is in space domain  

    try
        outputSet = cell(nScramble, 1);
        for i = 1:nScramble
        
            randImg = round(255*rand(size(inputImgMagnitude)));
            randPhase = angle(fft2(double(randImg)));
    
            cmplxIm = inputImgMagnitude.*exp(1i.*randPhase);
            outputSet{i} = ifft2(cmplxIm, 'symmetric');
        end
    catch err
        disp('scrambleSet:Error  ');
        disp(err.message);
        disp(err.stack(1));
        disp(err.stack(2));
    end
end
