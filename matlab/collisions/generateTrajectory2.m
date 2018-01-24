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
deltaThtN = ( mod(Target.tht - geom.tht1 - thteps, 2*pi) - ...
              mod(startN.tht - geom.tht1 - thteps, 2*pi));

deltaPhi = (Target.phi - startP.phi); 

if ~isempty(geom.F1Pm) % if there is a motor map...
    
    % calculate the number of steps for each motor of each positioner
    % pull in the maps.  reverse move maps are flipped so that motion
    % in time always increases the index in Map.x
    Map.thtP =        geom.F1Pm;
    Map.thtN = fliplr(geom.F1Nm);
    Map.phiP =        geom.F2Pm;
    Map.phiN = fliplr(geom.F2Nm); % unnecessary in this context.
    
    n1bins = size(geom.F1Pm,2);
    n2bins = size(geom.F2Pm,2);

    %% start bins for P and N moves, unflipped maps.
    strtBin.thtP = (mod(startP.tht - geom.tht0 + thteps, 2*pi) - thteps - geom.map_range.tht(1))/geom.binWidth;
    strtBin.thtN = (geom.map_range.tht(2) - (mod(startN.tht - geom.tht0, 2*pi) + ADD2Pi))/geom.binWidth;
    strtBin.phiP  = (startP.phi - geom.map_range.phi(1))/geom.binWidth;
    strtBin.phiN  = n2bins - strtBin.phiP; % start positions are physically identical
    
    %% starting index (integer) is ceiling of starting bin (real) with mininmum value of 1.
    strtINDX.thtP = max(ceil(strtBin.thtP), 1);
    strtINDX.thtN = max(ceil(strtBin.thtN), 1);
    strtINDX.phiP = max(ceil(strtBin.phiP), 1);
    strtINDX.phiN = max(ceil(strtBin.phiN), 1);
    
    % here, thteps takes care of near zero tht targets while leaving those too close to the opp-sense HS
    % for the same-sense situation.
    TGT_IN_OVERLAP = mod(Target.tht - geom.tht0 + thteps, 2*pi) < geom.tht_overlap;
    
    fnshBin.thtP = (mod(Target.tht - geom.tht0 + thteps, 2*pi) - thteps - geom.map_range.tht(1))/geom.binWidth;
    fnshBin.thtN = n1bins - fnshBin.thtP - 2*pi/geom.binWidth * TGT_IN_OVERLAP;
    fnshBin.phiP = (Target.phi - geom.map_range.phi(1))/geom.binWidth;
    fnshBin.phiN = n2bins - fnshBin.phiP;
% $$$     %% fix for roundoff errors near zero
% $$$     fnshBin.thtP = max(fnshBin.thtP,thteps);
% $$$     fnshBin.thtN = max(fnshBin.thtN,thteps);
% $$$     fnshBin.phiP = max(fnshBin.phiP,thteps);
% $$$     fnshBin.phiN = max(fnshBin.phiN,thteps);
    
    fnshINDX.thtP = max(ceil(fnshBin.thtP), 1);
    fnshINDX.thtN = max(ceil(fnshBin.thtN), 1);
    fnshINDX.phiP = max(ceil(fnshBin.phiP), 1);
    fnshINDX.phiN = max(ceil(fnshBin.phiN), 1);

    %% don't let the fnshINDX exceed n1bins
    tempbool = fnshINDX.thtN > n1bins;
    if sum(tempbool) > 0
        if sum(fnshBin.thtP(tempbool) < -thteps)
            disp(['Warning: fnshINDX exceeds n1bins in some cases probably due to roundoff.  ' ...
                  'Correcting positioners (don''t worry if thtp is small...'])
            disp(sprintf('pid %02d, thtp = %g\n',[find(tempbool) fnshBin.thtP(tempbool)].'))
        end
        fnshINDX.thtN(tempbool) = n1bins;
    end
        
    %% calculate the fractional bins at the Randbedingung of the trajectory
    strtOverCount.thtP = strtBin.thtP - (strtINDX.thtP - 1);
    strtOverCount.thtN = strtBin.thtN - (strtINDX.thtN - 1);
    strtOverCount.phiP = strtBin.phiP - (strtINDX.phiP - 1);
    strtOverCount.phiN = strtBin.phiN - (strtINDX.phiN - 1);
    fnshOverCount.thtP = fnshINDX.thtP - fnshBin.thtP;
    fnshOverCount.thtN = fnshINDX.thtN - fnshBin.thtN;
    fnshOverCount.phiP = fnshINDX.phiP - fnshBin.phiP;
    fnshOverCount.phiN = fnshINDX.phiN - fnshBin.phiN;

    % calculate number of steps required for each positioner
    for jj = 1:nCobras
        try
        nSteps.thtP(1,jj) = (sum(Map.thtP(jj, strtINDX.thtP(jj):fnshINDX.thtP(jj))) ...
                             - Map.thtP(jj,strtINDX.thtP(jj)) * strtOverCount.thtP(jj) ...
                             - Map.thtP(jj,fnshINDX.thtP(jj)) * fnshOverCount.thtP(jj) );
        nSteps.thtN(1,jj) = (sum(Map.thtN(jj, strtINDX.thtN(jj):fnshINDX.thtN(jj))) ...
                             - Map.thtN(jj,strtINDX.thtN(jj)) * strtOverCount.thtN(jj) ...
                             - Map.thtN(jj,fnshINDX.thtN(jj)) * fnshOverCount.thtN(jj) );
        nSteps.phiP(1,jj) = (sum(Map.phiP(jj, strtINDX.phiP(jj):fnshINDX.phiP(jj))) ...
                             - Map.phiP(jj,strtINDX.phiP(jj)) * strtOverCount.phiP(jj) ...
                             - Map.phiP(jj,fnshINDX.phiP(jj)) * fnshOverCount.phiP(jj) );
        nSteps.phiN(1,jj) = (sum(Map.phiN(jj, strtINDX.phiN(jj):fnshINDX.phiN(jj))) ...
                             - Map.phiN(jj,strtINDX.phiN(jj)) * strtOverCount.phiN(jj) ...
                             - Map.phiN(jj,fnshINDX.phiN(jj)) * fnshOverCount.phiN(jj) );
        catch
            disp(['generateTrajectory2 warning: Usually, when things ' ...
                  'go wrong, it is because the XML file is bad. Error in pid ' ...
                  num2str(geom.pids(jj))]);
        end
    end
    
    nSteps.thtP = max(nSteps.thtP, 0);
    nSteps.thtN = max(nSteps.thtN, 0);
    nSteps.phiP = max(nSteps.phiP, 0);
    nSteps.phiN = max(nSteps.phiN, 0);

    %% need to write Tht and Phi
    
    %% fractionalBinError derived in PHM's COO notebook #2, 1-dec-2015
    fractionalBinError = geom.alpha * geom.binWidth^(geom.beta - 1);

    %% Map has dimensions (#fibers,steps) where steps = 100 (tht),
    %% 50 (phi)

    % need to calculate a map error factor.  fBE * randn = DX must be
    % repeatedly applied until the cumulative sum of 1+DX exceeds 1
    % (the normalized bin width).  The error factor multiplied into
    % Map to make noisyMap is ((#-1) of random numbers + fraction of
    % the last one)
    noisyMap.thtP = Map.thtP .* mapFactor(fractionalBinError,size(Map.thtP));
    noisyMap.thtN = Map.thtN .* mapFactor(fractionalBinError,size(Map.thtN));
    noisyMap.phiP = Map.phiP .* mapFactor(fractionalBinError,size(Map.phiP));
    noisyMap.phiN = Map.phiN .* mapFactor(fractionalBinError,size(Map.phiN));
    
    % add an sticky bin at the end to prevent over-runs.
    noisyMap.thtP = [noisyMap.thtP ones(nCobras,1)*1e9];
    noisyMap.thtN = [noisyMap.thtN ones(nCobras,1)*1e9];
    noisyMap.phiP = [noisyMap.phiP ones(nCobras,1)*1e9];
    noisyMap.phiN = [noisyMap.phiN ones(nCobras,1)*1e9];

    stepsPerTime = 50;

    for jj = 1:nCobras
% $$$         %% It's a bad idea to assume we start from the first bin in the cumulative sum.
% $$$         trajEndBin = find( nSteps.thtP(jj) < cumsum(noisyMap.thtP(jj,:)), 1);
% $$$         stepCtr = [0, cumsum(noisyMap.thtP(jj, 1:trajEndBin))];
% $$$         binCtr  = 0:(length(stepCtr)-1);
% $$$         trajSteps  = [0:stepsPerTime:nSteps.thtP(jj), nSteps.thtP(jj)];
% $$$         if trajSteps(1) < stepCtr(1) | trajSteps(end) > stepCtr(end)
% $$$             fprintf(1,'Warning (genTraj2) extrapolating to determine +tht trajectory on fiber %d\n',jj);
% $$$             keyboard
% $$$         end
% $$$         ThtP{jj} = startP.tht(jj) + geom.binWidth * interp1(stepCtr, binCtr, trajSteps,'linear','extrap');
% $$$         tBins.thtP(jj) = length(trajSteps);

        %%% THT P (tht moving out)
        %% the motion is NOT assumed to start at angle = zero, bin = 1
        stepOffset = strtOverCount.thtP(jj) * noisyMap.thtP(jj,strtINDX.thtP(jj));
        trajEndBin = (find(nSteps.thtP(jj) < cumsum(noisyMap.thtP(jj,strtINDX.thtP(jj):end)) - stepOffset, 1)...
                      + strtINDX.thtP(jj) - 1);
        stepCtr = [0, cumsum(noisyMap.thtP(jj, strtINDX.thtP(jj):trajEndBin)) - stepOffset];
        %% there may be Randbedingung issues for non-moves.
        binCtr  = [strtBin.thtP(jj) strtINDX.thtP(jj) + (0:(length(stepCtr) - 2))];%0:(length(stepCtr)-1);
        trajSteps  = [0:stepsPerTime:nSteps.thtP(jj), nSteps.thtP(jj)];
        if trajSteps(1) < stepCtr(1) | trajSteps(end) > stepCtr(end)
            fprintf(1,'Warning (genTraj2) extrapolating to determine tht trajectory on fiber %d\n',jj);
        end
        ThtP{jj} = startP.tht(jj) + geom.binWidth * interp1(stepCtr, binCtr, trajSteps,'linear','extrap');
        tBins.thtP(jj) = length(trajSteps);

        %%% THT N (opp sense)
        stepOffset = strtOverCount.thtN(jj) * noisyMap.thtN(jj,strtINDX.thtN(jj));
        trajEndBin = (find(nSteps.thtN(jj) < cumsum(noisyMap.thtN(jj,strtINDX.thtN(jj):end)) - stepOffset, 1)...
                      + strtINDX.thtN(jj) - 1);
        stepCtr = [0, cumsum(noisyMap.thtN(jj, strtINDX.thtN(jj):trajEndBin)) - stepOffset];
        binCtr  = [strtBin.thtN(jj) strtINDX.thtN(jj) + (0:(length(stepCtr) - 2))];
        trajSteps  = [0:stepsPerTime:nSteps.thtN(jj), nSteps.thtN(jj)];
        if trajSteps(1) < stepCtr(1) | trajSteps(end) > stepCtr(end)
            fprintf(1,'Warning (genTraj2) extrapolating to determine -tht trajectory on fiber %d\n',jj);
        end
% $$$         ThtN{jj} = geom.tht0(jj) + geom.binWidth * (n1bins - interp1(stepCtr, binCtr, trajSteps,'linear','extrap'));
        ThtN{jj} = geom.tht0(jj) + geom.map_range.tht(end) - geom.binWidth * interp1(stepCtr, binCtr, trajSteps,'linear','extrap');
        tBins.thtN(jj) = length(trajSteps);

        %%% PHI P (phi moving out)
        stepOffset = strtOverCount.phiP(jj) * noisyMap.phiP(jj,strtINDX.phiP(jj));
        trajEndBin = (find(nSteps.phiP(jj) < cumsum(noisyMap.phiP(jj,strtINDX.phiP(jj):end)) - stepOffset, 1)...
                      + strtINDX.phiP(jj) - 1);
        stepCtr = [0, cumsum(noisyMap.phiP(jj, strtINDX.phiP(jj):trajEndBin)) - stepOffset];
        binCtr  = [strtBin.phiP(jj) strtINDX.phiP(jj) + (0:(length(stepCtr) - 2))];%0:(length(stepCtr)-1);
        trajSteps  = [0:stepsPerTime:nSteps.phiP(jj), nSteps.phiP(jj)];
        if trajSteps(1) < stepCtr(1) | trajSteps(end) > stepCtr(end)
            fprintf(1,'Warning (genTraj2) extrapolating to determine phi trajectory on fiber %d\n',jj);
        end
        PhiP{jj} = geom.map_range.phi(1) + geom.binWidth * interp1(stepCtr, binCtr, trajSteps,'linear','extrap');
        tBins.phiP(jj) = length(trajSteps);

        %%% PHI N (opp sense)  :::: NOT QUALIFIED
        stepOffset = strtOverCount.phiN(jj) * noisyMap.phiN(jj,strtINDX.phiN(jj));
        trajEndBin = (find(nSteps.phiN(jj) < cumsum(noisyMap.phiN(jj,strtINDX.phiN(jj):end)) - stepOffset, 1)...
                      + strtINDX.phiN(jj) - 1);
        stepCtr = [0, cumsum(noisyMap.phiN(jj, strtINDX.phiN(jj):trajEndBin)) - stepOffset];
        binCtr  = [strtBin.phiN(jj) strtINDX.phiN(jj) + (0:(length(stepCtr) - 2))];
        trajSteps  = [0:stepsPerTime:nSteps.phiN(jj), nSteps.phiN(jj)];
        if trajSteps(1) < stepCtr(1) | trajSteps(end) > stepCtr(end)
            fprintf(1,'Warning (genTraj2) extrapolating to determine -phi trajectory on fiber %d\n',jj);
        end
        PhiN{jj} = geom.map_range.phi(end) - geom.binWidth * interp1(stepCtr, binCtr, trajSteps,'linear','extrap');
        tBins.phiN(jj) = length(trajSteps);

    end

    %%% Positive and Negative direction trajectories now defined for both Tht and Phi stages
    %%% In valid trajectories, the index (trajSteps?) is monotonically increasing.  Sort out the valid Phi
    %%% trajectories and combine PhiP and PhiN so that all trajectories are valid.  Maintain book-keeping
    %%% on all of the helper variables (# steps and # bins)
    
    phi_moves_N = fnshBin.phiP < strtBin.phiP;

    Phi = PhiP;
    Phi(phi_moves_N) = PhiN(phi_moves_N);
    nSteps.phi = nSteps.phiP;
    nSteps.phi(phi_moves_N) = nSteps.phiN(phi_moves_N);
    tBins.phi = tBins.phiP;
    tBins.phi(phi_moves_N) = tBins.phiN(phi_moves_N);
    
    % should this be max(min(thtp, thtn), phi) ? see note below. 
    % nSteps.max is not presently an important variable.
    nSteps.max = max([nSteps.thtP; nSteps.thtN; nSteps.phi]);
    %%%
    
    % this is the longest vector needed for any trajectory.
    tBins.max = max([tBins.thtP tBins.thtN tBins.phi]);
end % there is always a motor map now    

output = packstruct(ThtP, ThtN, Phi);
output.nthtP = nSteps.thtP;
output.nthtN = nSteps.thtN;
output.nphi  = nSteps.phi;
output.nmax  = nSteps.max; % presently (1/24/18) not used
                           % downstream except for its size.
output.lthtP = tBins.thtP;
output.lthtN = tBins.thtN;
output.lphi  = tBins.phi;
output.lmax  = tBins.max;

%% verification plots
if exist('verify','var')
    
% $$$     disp(['[dphi dtht]: difference between target and final ' ...
% $$$         'trajectory position']);
% $$$     [mod(Tht(:,end) - Targets.tht + pi,2*pi) - pi, Phi(:,end) - Targets.phi]
    
    for jj = 1:length(output.nmax)
       thtP(jj,:) = padarray(ThtP{jj}, [0, tBins.max - length(ThtP{jj})], 'replicate','post');
       thtN(jj,:) = padarray(ThtN{jj}, [0, tBins.max - length(ThtN{jj})], 'replicate','post');
       phiP(jj,:) = padarray(PhiP{jj}, [0, tBins.max - length(PhiP{jj})], 'replicate','pre');
       phiN(jj,:) = padarray(PhiN{jj}, [0, tBins.max - length(PhiN{jj})], 'replicate','pre');
    end
    
    TrajP = (bsxfun(@times, geom.L1, exp(1i*thtP)) + ...
            bsxfun(@times, geom.L2, exp(1i*(thtP+phiP)) ) );
    TrajP = bsxfun(@plus, geom.center, TrajP);
    TrajN = (bsxfun(@times, geom.L1, exp(1i*thtN)) + ...
            bsxfun(@times, geom.L2, exp(1i*(thtN+phiP)) ) );
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
