function output=generateTrajectory(currentPosition,targetList, geom, trajectory_strategy, verify)
% generate a trajectory matrix given an assigned target list and
% center positions

if ~exist('trajectory_strategy','var'), trajectory_strategy = 'latelate'; end;

stepSize = 100e-3; %radians

% theta/phi positions of the current and target positions
strtPos = XY2TP(currentPosition - geom.center, geom.L1, geom.L2);
Targets = XY2TP(targetList - geom.center, geom.L1, geom.L2);

% some useful variables on the way to calculating the number of
% steps required to complete the trajectory.
deltaTht = mod(Targets.tht - strtPos.tht, 2*pi); % this is > 0 for fwd
deltaPhi = (Targets.phi - strtPos.phi); % this is always < 0 for fwd

thtFWD = deltaTht > 0;
phiFWD = deltaPhi < 0;

if isfield(geom,'S1Fm')
    
    % calculate the number of steps for each motor of each positioner
    % pull in the maps.  reverse move maps are flipped so that motion
    % in time always increases the index in Map.x
    n1bins = 100;
    n2bins = 50;;
    
    Map.tht = fliplr(geom.S1Rm(:,1:n1bins));
    Map.phi = fliplr(geom.S2Rm(:,1:n2bins));
    Map.tht(thtFWD,:) = geom.S1Fm(thtFWD,1:n1bins);
    Map.phi(phiFWD,:) = geom.S2Fm(phiFWD,1:n2bins);
    
    % real number bin indices - reverse movement indices handled by second term
    strtBin.tht = abs((mod(strtPos.tht-geom.tht0 + 1e-9, 2*pi)-1e-9)/geom.binWidth - ~thtFWD*n1bins);
    fnshBin.tht = abs((mod(Targets.tht-geom.tht0 + 1e-9, 2*pi)-1e-9)/geom.binWidth - ~thtFWD*n1bins);
    strtBin.phi = abs(strtPos.phi/geom.binWidth - phiFWD*n2bins);
    fnshBin.phi = abs(Targets.phi/geom.binWidth - phiFWD*n2bins);
    
    % calculate number of steps required for each positioner
    for jj = 1:length(geom.center)
        jjstrt = ceil(strtBin.tht(jj));
        jjfnsh = ceil(fnshBin.tht(jj));
        % start with total including the partial bins and then trim the ends
        nSteps.tht(jj) = sum(Map.tht(jj, jjstrt:jjfnsh)) ...
            - Map.tht(jj, jjstrt) * mod(strtBin.tht(jj),1) ...
            - Map.tht(jj, jjfnsh) * mod(-fnshBin.tht(jj),1) ;
        if nSteps.tht(jj) < 0
            disp('nsteps.tht < 0 detected');
            keyboard;
        end
        
        jjstrt = ceil(strtBin.phi(jj));
        jjfnsh = ceil(fnshBin.phi(jj));
        % start with total including the partial bins and then trim the ends
        nSteps.phi(jj) = sum(Map.phi(jj, jjstrt:jjfnsh)) ...
            - Map.phi(jj, jjstrt) * mod(strtBin.phi(jj),1) ...
            - Map.phi(jj, jjfnsh) * mod(-fnshBin.phi(jj),1);
        
    end
    nSteps.max = max([nSteps.tht nSteps.phi]);
    
    %% need to write Tht and Phi
    
    alpha = 0.07; % error of a 1 radian move
    beta  = 0.50; % scaling factor for move size
    
    fractionalBinError = alpha / geom.binWidth^beta;
    
    noisyMap.tht = Map.tht .* (1./ abs(1 + fractionalBinError * randn(size(Map.tht))));
    noisyMap.phi = Map.phi .* (1./ abs(1 + fractionalBinError * randn(size(Map.phi))));
    
    stepsPerTime = 50;
    lengthTime = ceil(nSteps.max/stepsPerTime);
    
    switch trajectory_strategy
      case 'earlyLate'
        thtPad = 'post';
        phiPad = 'pre';
      case 'lateLate'
        thtPad = 'pre';
        phiPad = 'pre';
      case 'earlyEarly'
        thtPad = 'post';
        phiPad = 'post';
    end
    
    
    for jj = 1:length(geom.center)
        jjstrt = ceil(strtBin.tht(jj));
        jjfnsh = ceil(fnshBin.tht(jj));
        frntExtraSteps = noisyMap.tht(jj,jjstrt) * mod(strtBin.tht(jj),1);
        trajEndBin = find( (nSteps.tht(jj) + frntExtraSteps) < cumsum(noisyMap.tht(jj,jjstrt:end)),1)...
            + jjstrt - 1;
        % backExtraSteps = sum(noisyMap.tht(jj,jjstrt:trajEndBin)) - frntExtraSteps - nSteps.tht(jj);
        
        stepCtr = [0, cumsum(noisyMap.tht(jj, jjstrt:trajEndBin))];
        binCtr  = 0:(length(stepCtr)-1);
        trajSteps  = [0:stepsPerTime:nSteps.tht(jj), nSteps.tht(jj)] + frntExtraSteps;
        Tht(jj,:) = padarray(interp1(stepCtr, binCtr, trajSteps) * geom.binWidth + strtPos.tht(jj), ...
                             [0 lengthTime - length(trajSteps) + 1], 'replicate',thtPad);

        %% theta has been checked out, but only for forward moves.  Phi is untested as of now.
        jjstrt = ceil(strtBin.phi(jj));
        jjfnsh = ceil(fnshBin.phi(jj));
        frntExtraSteps = noisyMap.phi(jj,jjstrt) * mod(strtBin.phi(jj),1);
        trajEndBin = find( (nSteps.phi(jj) + frntExtraSteps) < cumsum(noisyMap.phi(jj,jjstrt:end)),1)...
            + jjstrt - 1;
        % backExtraSteps = sum(noisyMap.phi(jj,jjstrt:trajEndBin)) - frntExtraSteps - nSteps.phi(jj);
        
        stepCtr = [0, cumsum(noisyMap.phi(jj, jjstrt:trajEndBin))];
        binCtr  = 0:(length(stepCtr)-1);
        trajSteps  = [0:stepsPerTime:nSteps.phi(jj), nSteps.phi(jj)] + frntExtraSteps;
        Phi(jj,:) = padarray(strtPos.phi(jj) - interp1(stepCtr, binCtr, trajSteps) * geom.binWidth, ...
                             [0 lengthTime - length(trajSteps) + 1], 'replicate',phiPad);
    end

% $$$     figure(1)
% $$$     plot(Tht'/(2*pi)); xlabel('time [50steps]'); ylabel('\Theta/(2*pi) [rad]');
% $$$     figure(2)
% $$$     plot(Phi'/pi); xlabel('time [50steps]'); ylabel('\Phi/\pi [rad]');
    
%     Plot for johannes to check credibility of noisy map.(something is still fishy with the values but i am off to lunch now).
     for kk = 1:length(geom.center)
         q = noisyMap.tht(kk,:);
         figure()
         hold on;
         qm = repmat(q,length(q),1);
         sumsq= sum(tril(qm),2);
         movstart = linspace(0,0,100)';
         movend = linspace(0,100,100)';
         movvel = movend./sumsq *180/pi;
         plot(movstart, movvel, 'rx');
         plot(movend, movvel, 'rx')
         plot([movstart, movend]', [movvel, movvel]', 'r')
     end
    keyboard;
    %%
else
    maxDeltaAngle = max(max(deltaTht), max(-deltaPhi));
    nSteps.max = ceil(maxDeltaAngle/stepSize);
    
    %% prototype for trajectories
    protoTraj = (0:nSteps.max) * stepSize;
    %  protoTraj = 0:stepSize:maxDeltaAngle;
    dmaxDeltaAngle = protoTraj(end) - maxDeltaAngle;
    
    
    %% define trajectories in theta/phi space
    Tht = bsxfun(@min, protoTraj, deltaTht);
    Tht = bsxfun(@plus, Tht, geom.tht0);
    
    switch trajectory_strategy
        case 'earlyLate'
            % the phi trajectory is set up so that it moves as late as possible
            Phi = bsxfun(@min , protoTraj, abs(deltaPhi));
            Phi = bsxfun(@plus, Phi      , Targets.phi);
            Phi = fliplr(Phi);
        case 'lateLate'
            % all movements end at on final step
            Tht = bsxfun(@plus, protoTraj - maxDeltaAngle - dmaxDeltaAngle, deltaTht);
            Tht = max(Tht,0);
            Tht = bsxfun(@plus, Tht, geom.tht0);
            Phi = bsxfun(@minus, protoTraj - maxDeltaAngle - dmaxDeltaAngle, deltaPhi);
            Phi = max(Phi,0);
            Phi = bsxfun(@minus, geom.phiIn, Phi);
        case 'earlyEarly' % no good!  don't use.  only here for historical reasons
            % phi moves right away
            Phi = bsxfun(@max, -protoTraj, deltaPhi);
            Phi = bsxfun(@plus, Phi, geom.phiIn);
    end
end

%% translate back into bench XY coordinates
Trajectories = (bsxfun(@times, geom.L1, exp(1i*Tht)) + ...
    bsxfun(@times, geom.L2, exp(1i*(Tht+Phi)) ) );
Trajectories = bsxfun(@plus, geom.center, Trajectories);

output = Trajectories;


%% verification plots
if exist('verify','var')
    
    disp(['[dphi dtht]: difference between target and final ' ...
        'trajectory position']);
    [mod(Tht(:,end) - Targets.tht + pi,2*pi) - pi, Phi(:,end) - Targets.phi]
    
    plot(Trajectories.');
    cmplx(@plotcircle,geom.center, geom.L1+geom.L2, 'k:');
    hold on;
    plot(targetList,'ro','MarkerFace','r');
    plot(Trajectories(:,end),'kx');
    plot(geom.L1.*exp(i*geom.tht0) + geom.L2.*exp(i*(geom.tht0+geom.phiIn)) + ...
        geom.center,'go','MarkerFace','g');
    
    
end