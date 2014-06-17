function avgMag = averageMagnitude(inputSet, domain)
% calculates average magnitude for input set
% magnitude will be calculated and returned in FREQUENCY domain  

    numEntries = length(inputSet);    
    avgMag = 0;
    try
        switch domain
            case 'FREQ'
                for i = 1:numEntries
                    avgMag = avgMag + abs(inputSet{i});  
                end
            case 'SPACE'
                for i = 1:numEntries
                    avgMag = avgMag + abs(fft2(inputSet{i}));  
                end
            otherwise
                warning('Unexpected domain, returning zero');
        end
         
    %the only possible error is dimentions mismatch
    catch err
        disp('averageMagnitude:Error  ');
        disp(err.message);
        disp(err.stack(1));
        disp(err.stack(2));
    end        
end