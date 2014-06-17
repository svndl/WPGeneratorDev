function scrambleSet(inputSet)
    domain = 'SPACE';
    
    numGroups = length(inputSet);
    result.s = cell(numGroups, 1);
    result.ps = cell(numGroups, 1);
    result.psA = cell(numGroups, 1);
    
    for i = 1:numGroups
        group = inputSet{i};
        [result.s{i}, result.ps{i}, result.psA{i}] = ScrambleBundle(group, domain);
    end
    save(strcat('result', datestr(clock), '.mat'), 'result');
end