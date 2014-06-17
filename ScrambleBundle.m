function [psScrambled, psScrambledAveraged, phaseScrambled] = ScrambleBundle(inputSet, domain)
    % scramble set in FREQ/SPACE domain
    
    numElements = length(inputSet);
    psScrambled = cell(numElements, 1);
    psScrambledAveraged = cell(numElements, 1);
    phaseScrambled = cell(numElements, 1);
    alteredInputSet = cell(numElements, 1);
    avgMag = averageMagnitude(inputSet, domain);
    numScrambled = 5;
    
    try
        switch domain
            case 'SPACE'
                for i = 1:numElements
                    [psScrambled{i}, alteredInputSet{i}] = psScramble(inputSet{i});
                    phaseScrambled{i} = scramblePhase(avgMag, numScrambled);
                end
            case 'FREQ'
                for i = 1:numElements
                    [psScrambled{i}, alteredInputSet{i}] = psScramble(ifft2(inputSet{i}));
                    phaseScrambled{i} = scramblePhase(avgMag, numScrambled);
                end
            otherwise
                disp('Domain is not specified, attempting FREQ');
                [psScrambled, psScrambledAveraged, phaseScrambled] = ScrambleBundle(inputSet, 'FREQ');
        end;
        alteredAvgMag = averageMagnitude(alteredInputSet, domain);
        newMag = avgMag;
        if (sum(size(alteredAvgMag) ~= size(avgMag))<2)
            newMag = alteredAvgMag;
        end;
        psScrambledAveraged = replaceMagnitude(psScrambled, newMag, 'SPACE');
    catch err
        disp('ScrambleBundle:Error  ');
        disp(err.message);
        disp(err.stack(1));
        disp(err.stack(2));
        psScrambled = 0;
        psScrambledAveraged = 0;
        phaseScrambled = 0;
    end
end

