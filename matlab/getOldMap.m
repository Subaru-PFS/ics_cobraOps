function output = getOldMap(stage, direction, id, CobraConfig)
% stage = 1/2,  direction  : 'fwd' or 'rev'
    fwdMapS1 = [];
    [speed angles] = getMMap(CobraConfig, id,'SLOW', stage, direction);
    fwdMapS1(1,:) = angles(3:end);
    fwdMapS1(2,:) = speed(3:end); 
    fwdMapS1(:,fwdMapS1(2,:)==0) = [];
    
    output = fwdMapS1;
end