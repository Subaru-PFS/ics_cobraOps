function noisy=addNoise(cumMap, fBE)
% add noise to a cumulative map, according to the fractionBinError (fBE).

    dataAxis = 2; % each cobras data is laid out in the 2nd axis

    Map = diff(cumMap,1,dataAxis);
    mapSize = size(Map);

    noisy = Map .* mapFactor(fBE, mapSize);
    noisy = cumsum([zeros(mapSize(1),1) Map], dataAxis);
    % make last entry huge to eliminate over-runs
    noisy(:,end) = 1e9;
end


%% eventually, mapFactor should go here.