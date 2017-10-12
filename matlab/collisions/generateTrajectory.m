function output=generateTrajectory(currentPosition,targetList, geom, trajectory_strategy, verify)
% generate a trajectory matrix given an assigned target list and
% center positions
%
% example: after running q = simFun(1.5,'full',1,1) try this:
% generateTrajectory(q.Traj.traj(:,end), q.targets, q.bench, 'earlyLate', 1);

if ~exist('trajectory_strategy','var'), trajectory_strategy = 'lateLate'; end;

stepSize = 100e-3; %radians, for non-motormap simulations

% this is an epsilon to make sure that values near zero in theta are
% intepreted on the positive (or negative) side of the cut,
% respectively. Make it negative for same same direction moveouts (positive
% for positive hardstop).
thteps = 1e-12;

nCobras = length(geom.center);

% theta/phi positions of the current and target positions
strtPos = XY2TP(currentPosition - geom.center, geom.L1, geom.L2);
Targets = XY2TP(targetList      - geom.center, geom.L1, geom.L2);

% some useful variables on the way to calculating the number of
% steps required to complete the trajectory.

%% need to know values for delta theta that do not travese a hard stop.
% $$$ % short and long are the directed paths from strtPos to Target
% $$$ deltatht.short = mod(Targets.tht - strtPos.tht + pi, 2*pi) - pi; % [-pi,pi]
% $$$ deltatht.long  = -(2*pi*sign(deltatht.short) - deltatht.short); % [-2pi,-pi],0,[pi,2pi]
% $$$ deltatht.from0 = mod(Targets.tht - geom.tht0, 2*pi) - mod(strtPos.tht - geom.tht0, 2*pi);
% $$$ deltatht.from1 = mod(Targets.tht - geom.tht1, 2*pi) - mod(strtPos.tht - geom.tht1, 2*pi);
% $$$ % It's tricky business to decide whether the long or the short path is the one that
% $$$ % doesn't cross a hard stop.  For now, let's assume that the short path is the one that
% $$$ % doesn't cross.  The only case I can think of right now is if the short path crosses BOTH
% $$$ % hard stops.  In that case, "from0" and "from1" will both match "long", so we'd have to
% $$$ % detect that condition.  In simulation, unless alpha is VERY large, that won't happen.

deltaTht = mod(Targets.tht - strtPos.tht + pi, 2*pi) - pi; % the "short" move

deltaPhi = (Targets.phi - strtPos.phi); 

active = abs(deltaTht) > thteps; % logical array, true for cobras that move.

% logical flag for forward trajectories
thtFWD = deltaTht > 0; 
phiFWD = deltaPhi > 0;
% direction is direction.
thtDIR = thtFWD * 2 - 1;
phiDIR = phiFWD * 2 - 1;

if ~isempty(geom.S1Pm) % if there is a motor map...
    
    % calculate the number of steps for each motor of each positioner
    % pull in the maps.  reverse move maps are flipped so that motion
    % in time always increases the index in Map.x
    n1bins = 100;
    n2bins = 50;
    
    Map.tht = fliplr(geom.S1Nm(:,1:n1bins));
    Map.phi = fliplr(geom.S2Nm(:,1:n2bins));
    Map.tht(thtFWD,:) = geom.S1Pm(thtFWD,1:n1bins);
    Map.phi(phiFWD,:) = geom.S2Pm(phiFWD,1:n2bins);
    
    % tht_stop is the hardstop position for the direction of movement
    tht_stop         = geom.tht1;
    tht_stop(thtFWD) = geom.tht0(thtFWD);
    
    % real number bin indices - reverse movement indices handled by second term
    strtBin.tht = thtDIR.*((mod(strtPos.tht-tht_stop + thteps, 2*pi)-thteps/2)/geom.binWidth - ~thtFWD*n1bins);
    fnshBin.tht = thtDIR.*((mod(Targets.tht-tht_stop + thteps, 2*pi)-thteps/2)/geom.binWidth - ~thtFWD*n1bins);
    try
        strtBin.phi = phiDIR.*((mod(strtPos.phi-geom.phiIn + pi/2, 2*pi) - pi/2)/geom.binWidth - ~phiFWD*n2bins);
        fnshBin.phi = phiDIR.*((mod(Targets.phi-geom.phiIn + pi/2, 2*pi) - pi/2)/geom.binWidth - ~phiFWD*n2bins);
    catch
        disp('generateTrajectory: args to mod were complex.  taking real part...')
        strtBin.phi = phiDIR.*((mod(real(strtPos.phi-geom.phiIn + pi/2), 2*pi) - pi/2)/geom.binWidth - ~phiFWD*n2bins);
        fnshBin.phi = phiDIR.*((mod(real(Targets.phi-geom.phiIn + pi/2), 2*pi) - pi/2)/geom.binWidth - ~phiFWD*n2bins);
    end

    % calculate number of steps required for each positioner
    for jj = 1:nCobras
        if active(jj)
            jjstrt = ceil(strtBin.tht(jj)); jjstrt = min(max(jjstrt,1),n1bins);
            jjfnsh = ceil(fnshBin.tht(jj)); jjfnsh = min(max(jjfnsh,1),n1bins);
            % start with total including the partial bins and then trim the ends
            try
                nSteps.tht(jj) = sum(Map.tht(jj, jjstrt:jjfnsh)) ...
                    - Map.tht(jj, jjstrt) * mod(strtBin.tht(jj),1) ...
                    - Map.tht(jj, jjfnsh) * mod(-fnshBin.tht(jj),1) ;
            catch
                keyboard
            end
            if nSteps.tht(jj) < 0
                disp('nsteps.tht < 0 detected');
            end
            % ceil because matlab starts at 1, not zero.
            jjstrt = ceil(strtBin.phi(jj)); jjstrt = min(max(jjstrt,1),n2bins);
            jjfnsh = ceil(fnshBin.phi(jj)); jjfnsh = min(max(jjfnsh,1),n2bins);
            % start with total including the partial bins and then trim the ends
            nSteps.phi(jj) = sum(Map.phi(jj, jjstrt:jjfnsh)) ...
                - Map.phi(jj, jjstrt) * mod(strtBin.phi(jj),1) ...
                - Map.phi(jj, jjfnsh) * mod(-fnshBin.phi(jj),1);
        else
            nSteps.tht(jj) = 0;
            nSteps.phi(jj) = 0;
        end
    end
    nSteps.max = max([nSteps.tht nSteps.phi]);
    nSteps.tht = nSteps.tht(:);
    nSteps.phi = nSteps.phi(:);

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
    noisyMap.tht = Map.tht .* mapFactor(fractionalBinError,size(Map.tht));
    noisyMap.phi = Map.phi .* mapFactor(fractionalBinError,size(Map.phi));
% $$$     noisyMap.tht = Map.tht .* (1./ abs(1 + fractionalBinError * randn(size(Map.tht))));
% $$$     noisyMap.phi = Map.phi .* (1./ abs(1 + fractionalBinError * randn(size(Map.phi))));
    
    % add an sticky bin at the end to prevent over-runs.
    noisyMap.tht = [noisyMap.tht ones(nCobras,1)*1e9];
    noisyMap.phi = [noisyMap.phi ones(nCobras,1)*1e9];

    stepsPerTime = 50;
    lengthTime = ceil(nSteps.max/stepsPerTime);
    
    switch trajectory_strategy
      case 'earlyLate'
        thtPad = 'post'; nSteps.dtht = zeros(nCobras,1);
        phiPad = 'pre';  nSteps.dphi = nSteps.max - nSteps.phi;
      case 'lateLate'
        thtPad = 'pre';  nSteps.dtht = nSteps.max - nSteps.tht;
        phiPad = 'pre';  nSteps.dphi = nSteps.max - nSteps.phi;
      case 'earlyEarly'
        thtPad = 'post'; nSteps.dtht = zeros(nCobras,1);
        phiPad = 'post'; nSteps.dphi = zeros(nCobras,1);       
      case 'lateEarly'
        thtPad = 'pre';  nSteps.dtht = nSteps.max - nSteps.tht;
        phiPad = 'post'; nSteps.dphi = zeros(nCobras,1);       
    end
    
    
    for jj = 1:nCobras
        jjstrt = ceil(strtBin.tht(jj));
        jjfnsh = ceil(fnshBin.tht(jj));
        try
            frntExtraSteps = noisyMap.tht(jj,jjstrt) * mod(strtBin.tht(jj),1);
        catch
            keyboard;
        end
        trajEndBin = find( (nSteps.tht(jj) + frntExtraSteps) < ...
                           cumsum(noisyMap.tht(jj,jjstrt:end)),1) + jjstrt - 1;
        stepCtr = [0, cumsum(noisyMap.tht(jj, jjstrt:trajEndBin))];
        binCtr  = 0:(length(stepCtr)-1);
        trajSteps  = [0:stepsPerTime:nSteps.tht(jj), nSteps.tht(jj)] + frntExtraSteps;
        if trajSteps(1) < stepCtr(1) | trajSteps(end) > stepCtr(end)
            fprintf(1,'Warning extrapolating to determine tht trajectory on fiber %d\n',jj);
        end
        Tht(jj,:) = padarray(strtPos.tht(jj) + thtDIR(jj) * geom.binWidth * ...
                             (interp1(stepCtr, binCtr, trajSteps,'linear','extrap') - ...
                              mod(strtBin.tht(jj),1)),...
                             [0 lengthTime - length(trajSteps) + 1], 'replicate',thtPad);
        %% theta has been checked out, but only for forward moves.  Phi is untested as of now.
        jjstrt = ceil(strtBin.phi(jj)); jjstrt = min(max(jjstrt,1),n2bins);
        jjfnsh = ceil(fnshBin.phi(jj)); jjfnsh = min(max(jjfnsh,1),n2bins);
        frntExtraSteps = noisyMap.phi(jj,jjstrt) * mod(strtBin.phi(jj),1);
        trajEndBin = find( (nSteps.phi(jj) + frntExtraSteps) < ...
                           cumsum(noisyMap.phi(jj,jjstrt:end)),1) + jjstrt - 1;
        stepCtr = [0, cumsum(noisyMap.phi(jj, jjstrt:trajEndBin))];
        binCtr  = 0:(length(stepCtr)-1);
        trajSteps  = [0:stepsPerTime:nSteps.phi(jj), nSteps.phi(jj)] + frntExtraSteps;
        % linear is default method for interp1.  'extrap' allows an
        % extrapolation outside of the bounds.  NOT a good idea,
        % except it helps handle some Phi situations.
        if trajSteps(1) < stepCtr(1) | trajSteps(end) > stepCtr(end)
            fprintf(1,'Warning extrapolating to determine phi trajectory on fiber %d\n',jj);
        end
        Phi(jj,:) = padarray(strtPos.phi(jj) + phiDIR(jj) * geom.binWidth * ...
                             (interp1(stepCtr, binCtr, trajSteps,'linear','extrap') -...
                              mod(strtBin.phi(jj),1)), ...
                             [0 lengthTime - length(trajSteps) + 1], 'replicate',phiPad);
        %        keyboard;
    end

% $$$     figure(1)
% $$$     plot(Tht'/(2*pi)); xlabel('time [50steps]'); ylabel('\Theta/(2*pi) [rad]');
% $$$     figure(2)
% $$$     plot(Phi'/pi); xlabel('time [50steps]'); ylabel('\Phi/\pi [rad]');
    
    %% Plot for johannes to check credibility of noisy map.(something is still fishy with the values but i am off to lunch now).
% $$$     for kk = 1:length(geom.center)
% $$$         q = noisyMap.tht(kk,:);
% $$$         figure()
% $$$         hold on;
% $$$         qm = repmat(q,length(q),1);
% $$$         sumsq= sum(tril(qm),2);
% $$$         movstart = linspace(0,0,100)';
% $$$         movend = linspace(0,100,100)';
% $$$         movvel = movend./sumsq;
% $$$         plot(movstart, movvel, 'rx');
% $$$         plot(movend, movvel, 'rx')
% $$$         plot([movstart, movend]', [movvel, movvel]', 'r')
% $$$     end
% $$$ keyboard;    
    %%
else % if there is no motor map...
    
    %%%%%  PHI REVERSE MOVES NOT IMPLEMENTED YET
    maxDeltaAngle = max(max(abs(deltaTht)), max(abs(deltaPhi)));
    nSteps.max = ceil(maxDeltaAngle/stepSize);
    nSteps.tht = ceil(deltaTht(:)/stepSize);
    nSteps.phi = ceil(deltaPhi(:)/stepSize);
    
    %% prototype for trajectories
    protoTraj = (0:nSteps.max) * stepSize;
    %  protoTraj = 0:stepSize:maxDeltaAngle;
    dmaxDeltaAngle = protoTraj(end) - maxDeltaAngle;
    
    
    %% define trajectories in theta/phi space
    %% early tht move defined here
    Tht = bsxfun(@min, protoTraj, abs(deltaTht));
    Tht = bsxfun(@times, Tht, thtDIR);
    Tht = bsxfun(@plus , Tht, strtPos.tht);
    %% late phi move
    Phi = bsxfun(@plus, protoTraj - maxDeltaAngle - dmaxDeltaAngle, abs(deltaPhi));
    Phi = max(Phi,0);
    Phi = bsxfun(@times, Phi, phiDIR);
    Phi = bsxfun(@plus,  Phi, strtPos.phi);
                                      
    switch trajectory_strategy                                      
        case 'earlyLate'
          nSteps.dtht = zeros(1,nCobras);       
          nSteps.dphi = nSteps.max - nSteps.phi;
        case 'lateLate'
            % all movements end at on final step
            Tht = bsxfun(@plus, protoTraj - maxDeltaAngle - dmaxDeltaAngle, abs(deltaTht));
            Tht = max(Tht,0);
            Tht = bsxfun(@times, Tht, thtDIR);
            Tht = bsxfun(@plus, Tht, geom.tht0);
            nSteps.dtht = nSteps.max - nSteps.tht;
            nSteps.dphi = nSteps.max - nSteps.phi;
        case 'earlyEarly' % no good!  don't use.  only here for historical reasons
            % phi moves right away
            Phi = bsxfun(@min, protoTraj, abs(deltaPhi));
            Phi = bsxfun(@times, Phi, phiDIR);
            Phi = bsxfun(@plus, Phi, strtPos.phi);
            nSteps.dtht = zeros(nCobras,1);       
            nSteps.dphi = zeros(nCobras,1);
    end
end

%% translate back into bench XY coordinates
Trajectories = (bsxfun(@times, geom.L1, exp(1i*Tht)) + ...
                bsxfun(@times, geom.L2, exp(1i*(Tht+Phi)) ) );
Trajectories = bsxfun(@plus, geom.center, Trajectories);

output.traj = Trajectories;
output.ntht = nSteps.tht;
output.dtht = nSteps.dtht;
output.nphi = nSteps.phi;
output.dphi = nSteps.dphi;
output.nmax = nSteps.max;

%% verification plots
if exist('verify','var')
    
% $$$     disp(['[dphi dtht]: difference between target and final ' ...
% $$$         'trajectory position']);
% $$$     [mod(Tht(:,end) - Targets.tht + pi,2*pi) - pi, Phi(:,end) - Targets.phi]
    
    tp = XY2TP(Trajectories - geom.center, geom.L1, geom.L2);
    elbow = geom.center + geom.L1.*exp(i*tp.tht);


    plot(Trajectories.','b');
    cmplx(@plotcircle,geom.center, geom.L1+geom.L2, 'k:');
    hold on;
    plot(elbow.','b:');
    pp(2) = plot(currentPosition,'go','MarkerFace','g');
    pp(3) = plot(Trajectories(:,end),'ro','MarkerFace','r');
    pp(1) = plot(targetList,'kx');
    plot([geom.center,  geom.center + (geom.L1+geom.L2).*exp(i*geom.tht0)].','k:');
    plot([geom.center,  geom.center + (geom.L1+geom.L2).*exp(i*geom.tht1)].','k:');
% $$$     plot(geom.L1.*exp(i*geom.tht0) + geom.L2.*exp(i*(geom.tht0+geom.phiIn)) + ...
% $$$         geom.center,'go','MarkerFace','g');
    legend(pp,'Target','Traj. start','Traj. end','Location','best');
end

return