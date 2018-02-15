function output=generateTrajectory(currentPosition,targetList, geom, verify)
% generate a trajectory matrix given an assigned target list and
% center positions
%
% example: after running q = simFun(1.5,'full',1,1) try this:
% generateTrajectory(q.Traj.traj(:,end), q.targets, q.bench, 1);

nCobras = length(geom.center);

% theta/phi positions of the current and target positions
strtPos = XY2TP(currentPosition - geom.center, geom.L1, geom.L2);
Targets = XY2TP(targetList      - geom.center, geom.L1, geom.L2);

strtTht0 = round(mod(strtPos.tht - geom.tht0, 2*pi), 13);
fnshTht0 = round(mod(Targets.tht - geom.tht0, 2*pi), 13);
strtTht1 = strtTht0 + 2*pi*(strtTht0 < geom.tht_overlap);
fnshTht1 = fnshTht0 + 2*pi*(fnshTht0 < geom.tht_overlap);

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

dtht0 = fnshTht0 - strtTht0;
dtht1 = fnshTht1 - strtTht1;
Use1  = (abs(dtht1) < abs(dtht0)); % strtTht or fnshTht is > 2*pi.
strtTht = strtTht0.*~Use1 + strtTht1.*Use1;
fnshTht = fnshTht0.*~Use1 + fnshTht1.*Use1;

deltaTht = fnshTht - strtTht; % the shortest available move that doens't cross boundaries.

deltaPhi = (Targets.phi - strtPos.phi); 

active = abs(deltaTht) > 0; % logical array, true for cobras that move.

% logical flag for forward trajectories
thtFWD = deltaTht >= 0; 
phiFWD = deltaPhi >= 0;

if ~isempty(geom.S1Pm) % if there is a motor map...
    
    % pull in the maps.
    Map.thtP = cumsum([zeros(nCobras,1) geom.F1Pm],2);
    Map.thtN = cumsum([zeros(nCobras,1) geom.F1Nm],2);
    Map.phiP = cumsum([zeros(nCobras,1) geom.F2Pm],2);
    Map.phiN = cumsum([zeros(nCobras,1) geom.F2Nm],2);
    
    n1bins = size(geom.F1Pm,2);
    n2bins = size(geom.F2Pm,2);

    ang.tht = (0:geom.binWidth:(geom.binWidth*n1bins))'; % tht indexes off SS hardstop
    ang.phi = (0:geom.binWidth:(geom.binWidth*n2bins))'-pi; %phi indexes off -pi

    %% # steps to command (interpolate strt/fnshTht against (ang,Map))
    for jj = 1:nCobras
        try
            if thtFWD(jj)
                nSteps.tht(1,jj) = (interp1(ang.tht, Map.thtP(jj,:), fnshTht(jj)) - ...
                                     interp1(ang.tht, Map.thtP(jj,:), strtTht(jj)));
            else
                nSteps.tht(1,jj) = (interp1(ang.tht, Map.thtN(jj,:), fnshTht(jj)) - ...
                                     interp1(ang.tht, Map.thtN(jj,:), strtTht(jj)));
            end
            if phiFWD(jj)
                nSteps.phi(1,jj) = (interp1(ang.phi, Map.phiP(jj,:), Targets.phi(jj)) - ...
                                     interp1(ang.phi, Map.phiP(jj,:), strtPos.phi(jj)));
            else
                nSteps.phi(1,jj) = (interp1(ang.phi, Map.phiN(jj,:), Targets.phi(jj)) - ...
                                     interp1(ang.phi, Map.phiN(jj,:), strtPos.phi(jj)));
            end
        catch
            disp(['generateTrajectory2 warning: Usually, when things ' ...
                  'go wrong, it is because the XML file is bad. Error in pid ' ...
                  num2str(geom.pids(jj))]);
        end
    end

    %% # noise up the map
    % fractionalBinError derived in PHM's COO notebook #2, 1-dec-2015
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

    %% find final position and trajectory with # steps - final location ~= target location.
    stepsPerTime = 50;
    for jj = 1:nCobras
        if thtFWD(jj)
            trajSteps  = ([0:stepsPerTime:nSteps.tht(jj), nSteps.tht(jj)] +...
                          interp1(ang.tht, noisyMap.thtP(jj,:), strtTht(jj)));
            protoTht{jj} = geom.tht0(jj) + interp1(noisyMap.thtP(jj,:), ang.tht, trajSteps);
        else
            trajSteps  = ([0:-stepsPerTime:nSteps.tht(jj), nSteps.tht(jj)] +...
                          interp1(ang.tht, noisyMap.thtN(jj,:), strtTht(jj)));
            protoTht{jj} = geom.tht0(jj) + interp1(noisyMap.thtN(jj,:), ang.tht, trajSteps);
        end
        tBins.tht(jj) = length(trajSteps);

        if phiFWD(jj)
            trajSteps  = ([0:stepsPerTime:nSteps.phi(jj), nSteps.phi(jj)] +...
                          interp1(ang.phi, noisyMap.phiP(jj,:), strtPos.phi(jj)));
            protoPhi{jj} = interp1(noisyMap.phiP(jj,:), ang.phi, trajSteps);
        else
            trajSteps  = ([0:-stepsPerTime:nSteps.phi(jj), nSteps.phi(jj)] +...
                          interp1(ang.phi, noisyMap.phiN(jj,:), strtPos.phi(jj)));
            protoPhi{jj} = interp1(noisyMap.phiN(jj,:), ang.phi, trajSteps);
        end
        tBins.phi(jj) = length(trajSteps);

    end
    % this is the longest vector needed for any trajectory.
    tBins.max = max([tBins.tht tBins.phi]);

    %% Turn Tht/Phi stubby arrays into uniform-length arrays.
    for jj = 1:nCobras
        if phiFWD(jj)
            thtPad = 'post'; % traditional earlyLate for outgoing trajectories.
            phiPad = 'pre';
        else
            thtPad = 'pre';  % run phi late if it's moving to the inside
            phiPad = 'post'; % run phi early if it's moving to the inside
        end
        Tht(jj,:) = padarray(protoTht{jj}, [0, tBins.max - tBins.tht(jj)],'replicate',thtPad);
        Phi(jj,:) = padarray(protoPhi{jj}, [0, tBins.max - tBins.phi(jj)],'replicate',phiPad);
    end
    
else % if there is no motor map...
    %% WARNING : as of 2/15/2018, this branch of the code is behind.

    stepSize = 100e-3; %radians, for non-motormap simulations

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

output = Trajectories;
% $$$ output.ntht = nSteps.tht;
% $$$ output.dtht = nSteps.dtht;
% $$$ output.nphi = nSteps.phi;
% $$$ output.dphi = nSteps.dphi;
% $$$ output.nmax = nSteps.max;

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