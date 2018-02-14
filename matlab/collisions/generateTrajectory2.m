function output=generateTrajectory2(targetList, geom, verify)
% Trajectory from home positions.
%
% "2" is for two hard stop options.

if ~exist('trajectory_strategy','var'), trajectory_strategy = 'lateLate'; end;

% this is an epsilon to make sure that values near zero in theta are
% intepreted on the positive (or negative) side of the cut,
% respectively. Make it negative for same same direction moveouts (positive
% for positive hardstop).
thteps = 1e-12; %geom.thteps;

nCobras = length(geom.center);

% theta/phi positions of the current and target positions
startP = XY2TP(geom.home0 - geom.center, geom.L1, geom.L2);
startN = XY2TP(geom.home1 - geom.center, geom.L1, geom.L2);
Target = XY2TP(targetList - geom.center, geom.L1, geom.L2);

% for relative angles, N direction moves need 2*pi added if they are in the overlap region
ADD2Pi = (mod(startN.tht - geom.tht0, 2*pi) < geom.tht_overlap + thteps) * 2 * pi;

% some useful variables on the way to calculating the number of
% steps required to complete the trajectory.
deltaThtP = ( mod(Target.tht - geom.tht0 + thteps, 2*pi) - ...
              mod(startP.tht - geom.tht0 + thteps, 2*pi));
deltaThtN = ( mod(startN.tht - geom.tht1 - thteps, 2*pi) - ...
              mod(Target.tht - geom.tht1 - thteps, 2*pi));

deltaPhi = (Target.phi - startP.phi); 

if ~isempty(geom.F1Pm) % if there is a motor map...
    
    %% calculate the number of steps for each motor of each positioner

    % pull in the maps.
    Map.thtP = cumsum([zeros(nCobras,1) geom.F1Pm],2);
    Map.thtN = cumsum([zeros(nCobras,1) geom.F1Nm],2);
    Map.phiP = cumsum([zeros(nCobras,1) geom.F2Pm],2);
    Map.phiN = cumsum([zeros(nCobras,1) geom.F2Nm],2);
    
    n1bins = size(geom.F1Pm,2);
    n2bins = size(geom.F2Pm,2);

    ang.tht = (0:geom.binWidth:(geom.binWidth*n1bins))'; % tht indexes off SS hardstop
    ang.phi = (0:geom.binWidth:(geom.binWidth*n2bins))'-pi; %phi indexes off -pi
    bins.tht = 0:n1bins;
    bins.phi = 0:n2bins;
    
    %% start bins for P and N moves, unflipped maps.
    strtBin.thtP = max(interp1(ang.tht,bins.tht, mod(startP.tht - geom.tht0 + thteps, 2*pi) - thteps),0);
    strtBin.thtN = min(interp1(ang.tht,bins.tht, mod(startN.tht - geom.tht0, 2*pi) + ADD2Pi), n1bins);
    strtBin.phi  = max((startP.phi - geom.map_range.phi(1))/geom.binWidth, 0);

    % here, thteps takes care of near zero tht targets while leaving those too close to the opp-sense HS
    % for the same-sense situation.
    TGT_IN_OVERLAP = mod(Target.tht - geom.tht0 + thteps, 2*pi) < geom.tht_overlap;
    
    fnshBin.thtP = interp1(ang.tht,bins.tht, mod(Target.tht - geom.tht0 - thteps, 2*pi) + thteps);
    fnshBin.thtN = interp1(ang.tht,bins.tht, mod(Target.tht - geom.tht0, 2*pi) + 2*pi*TGT_IN_OVERLAP);
    fnshBin.phi  = (Target.phi - geom.map_range.phi(1))/geom.binWidth;

    % calculate number of steps required for each positioner
    for jj = 1:nCobras
        try
            nSteps.thtP(1,jj) = (interp1(bins.tht, Map.thtP(jj,:), fnshBin.thtP(jj)) - ...
                                 interp1(bins.tht, Map.thtP(jj,:), strtBin.thtP(jj)));
            nSteps.thtN(1,jj) = (interp1(bins.tht, Map.thtN(jj,:), fnshBin.thtN(jj)) - ...
                                 interp1(bins.tht, Map.thtN(jj,:), strtBin.thtN(jj)));
            nSteps_phiP(1,jj) = (interp1(bins.phi, Map.phiP(jj,:), fnshBin.phi(jj)) - ...
                                 interp1(bins.phi, Map.phiP(jj,:), strtBin.phi(jj)));
            nSteps_phiN(1,jj) = (interp1(bins.phi, Map.phiN(jj,:), fnshBin.phi(jj)) - ...
                                 interp1(bins.phi, Map.phiN(jj,:), strtBin.phi(jj)));
        
        catch
            disp(['generateTrajectory2 warning: Usually, when things ' ...
                  'go wrong, it is because the XML file is bad. Error in pid ' ...
                  num2str(geom.pids(jj))]);
        end
    end

    nSteps.phi = max(nSteps_phiP, -nSteps_phiN);
    phiIsP = nSteps_phiP >= 0; % logical to track trajectories that move out in phi
    
    %% need to write Tht and Phi
    
    %% fractionalBinError derived in PHM's COO notebook #2, 1-dec-2015
    fractionalBinError = geom.alpha * geom.binWidth^(geom.beta - 1);

    % need to calculate a map error factor.  fBE * randn = DX must be
    % repeatedly applied until the cumulative sum of 1+DX exceeds 1
    % (the normalized bin width).  The error factor multiplied into
    % Map to make noisyMap is ((#-1) of random numbers + fraction of
    % the last one)

    % addNoise takes a cumsum Map and returns a cumsum Map with noise defined by fBE
    noisyMap.thtP = addNoise(Map.thtP, fractionalBinError);
    noisyMap.thtN = addNoise(Map.thtN, fractionalBinError);
    noisyMap.phiP = addNoise(Map.phiP, fractionalBinError);
    noisyMap.phiN = addNoise(Map.phiN, fractionalBinError);
    
    stepsPerTime = 50;

    for jj = 1:nCobras
        %%% THT P (tht moving out SS)
        trajSteps  = ([0:stepsPerTime:nSteps.thtP(jj), nSteps.thtP(jj)] +...
                      interp1(bins.tht, noisyMap.thtP(jj,:), strtBin.thtP(jj)));
        ThtP{jj} = geom.tht0(jj) + interp1(noisyMap.thtP(jj,:), ang.tht, trajSteps);
        tBins.thtP(jj) = length(trajSteps);
        %%% THT N (tht moving out OS)
        trajSteps  = ([0:-stepsPerTime:nSteps.thtN(jj), nSteps.thtN(jj)] +...
                      interp1(bins.tht, noisyMap.thtN(jj,:), strtBin.thtN(jj)));
        ThtN{jj} = geom.tht0(jj) + interp1(noisyMap.thtN(jj,:), ang.tht, trajSteps);
        tBins.thtN(jj) = length(trajSteps);

        if phiIsP(jj) % (nsteps.phiP(jj) >= 0)
            %%% PHI P (phi moving out)
            trajSteps  = ([0:stepsPerTime:nSteps.phi(jj), nSteps.phi(jj)] +...
                          interp1(bins.phi, noisyMap.phiP(jj,:), strtBin.phi(jj)));
            Phi{jj} = interp1(noisyMap.phiP(jj,:), ang.phi, trajSteps);
        else
            %%% PHI N (phi moving in)
            trajSteps  = ([0:-stepsPerTime:nSteps.phi(jj), nSteps.phi(jj)] +...
                          interp1(bins.phi, noisyMap.phiN(jj,:), strtBin.phi(jj)));
            Phi{jj} = interp1(noisyMap.phiN(jj,:), ang.phi, trajSteps);
        end
        tBins.phi(jj) = length(trajSteps);

    end

    %%% Positive and Negative direction trajectories now defined for both Tht and Phi stages
    
    % this is the longest vector needed for any trajectory.
    tBins.max = max([tBins.thtP tBins.thtN tBins.phi]);
end

output = packstruct(ThtP, ThtN, Phi);
output.nthtP = nSteps.thtP;
output.nthtN = nSteps.thtN;
output.nphi  = nSteps.phi;
output.lthtP = tBins.thtP;
output.lthtN = tBins.thtN;
output.lphi  = tBins.phi;
output.lmax  = tBins.max;

%% verification plots
if exist('verify','var')
    
% $$$     disp(['[dphi dtht]: difference between target and final ' ...
% $$$         'trajectory position']);
% $$$     [mod(Tht(:,end) - Targets.tht + pi,2*pi) - pi, Phi(:,end) - Targets.phi]
    
    for jj = 1:nCobras
       thtP(jj,:) = padarray(ThtP{jj}, [0, tBins.max - length(ThtP{jj})], 'replicate','post');
       thtN(jj,:) = padarray(ThtN{jj}, [0, tBins.max - length(ThtN{jj})], 'replicate','post');
       phi(jj,:) = padarray(Phi{jj}, [0, tBins.max - length(Phi{jj})], 'replicate','pre');
    end
    
    TrajP = (bsxfun(@times, geom.L1, exp(1i*thtP)) + ...
            bsxfun(@times, geom.L2, exp(1i*(thtP+phi)) ) );
    TrajP = bsxfun(@plus, geom.center, TrajP);
    TrajN = (bsxfun(@times, geom.L1, exp(1i*thtN)) + ...
            bsxfun(@times, geom.L2, exp(1i*(thtN+phi)) ) );
    TrajN = bsxfun(@plus, geom.center, TrajN);


    plot(TrajP.','b'); hold on;
    plot(TrajN.','c');
    cmplx(@plotcircle,geom.center, geom.L1+geom.L2, 'k:');
    plot(targetList,'ro','MarkerFace','r');
    plot(TrajN(:,end),'kx');
    plot(geom.L1.*exp(i*geom.tht0) + geom.L2.*exp(i*(geom.tht0+geom.phiIn)) + ...
        geom.center,'go','MarkerFace','g');
    hold off;
end

return
end
