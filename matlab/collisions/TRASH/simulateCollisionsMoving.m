% starter file for collision investigations

close all;
clear all;
showMoves = true; % Show single collision movements
%% Specify number of targets to generate
numtrg = 1000; 

%% Specify move strategy
%moveStrategy = 'earlyEarly';
%altmoveStrategy = 'earlyEarly';
moveStrategy = 'earlyLate'; % aka late, or L hook
altmoveStrategy = 'lateLate'; % aka radial, or R hook

useRealMaps = true;

%% Generate Positioner System
%% Q: why is the 7 pos system worse than the 14 rail system?

centers = [];

% 7 positioners in hex pattern      (R3 collisions = 3.0e-4 w/ 100k targets)
%centers = [0 8 * exp(1i * (0:5) * pi / 3)].';
 
% Line of cobras                    (R3 collisions = 4.5e-5 (fraction) w/ 1M targets)
centers = complex([0:8:(26*8)],0)' + i;

% N rails                           (R3 collisions = 9.9e-5 w/ 16k targets)
%centers = getCentersRails(14);

% full PFI bench                    (R3 collisions = 1.7e-4 (fraction) w/ 1M targets)
%centers = getCentersRails(14);
%centers = bsxfun(@times, [1 exp(i*2*pi/3) exp(i*4*pi/3)], centers);
%centers = reshape(centers,[],1);

%% Define a structure with the bench geometry

bench = defineBenchGeometry(centers,useRealMaps);

% reassign bench alpha/beta here (if desired)
bench.alpha=0.1;

numPos = length(bench.center);
numFields = ceil(numtrg / numPos);

%% Test new set of targets using rules:

% angular coordinate of target
THT = rand(numPos,numFields)*2*pi;

% radial coordinate of target
dA = 1 ./ ( (bench.rMax./bench.rMin).^2 - 1 ); %fraction of the keepout area.
rRange = sqrt(bench.rMax.^2 - bench.rMin.^2);
RDS = bsxfun(@times, sqrt(bsxfun(@plus, rand(size(THT)), dA)), rRange);

%% Rule 2 No targets closer than 2mm to any line in its target location. 
%disp('Number of bad targets from rule #2 (>2mm distance)')
targets   = bsxfun(@plus, RDS.*exp(1i*THT), bench.center);
targetsTP = XY2TP(bsxfun(@minus, targets, bench.center), bench.L1, bench.L2);

% $$$ figure(9999); plot(targets,'b.'); hold on; axis equal;
% $$$ title('targets'); xlabel('X [mm]'); ylabel('Y [mm]');

numcoll = 0;
collCenter = 0;
tic

for kk = 1:numFields
    firstRun = true;
    while(true)
        distances = CalcDistanceMatrix(targets(:,kk), bench);
        %        distances = symmetrize(distances);
        distances = min(distances, distances');
        [R2.row R2.col] = find(distances > 0 & distances < bench.minDist);
        replanCobra = R2.row(R2.row < R2.col); % only replan half
                                               % of the cobras.
% $$$         replanIndx = find(distances.dst < bench.minDist & ... % find colliders
% $$$                           distances.rc(:,1) < distances.rc(:,2)); %only replan the upper triangle
% $$$         replanCobra = distances.rc(replanIndx,1);
        nReplan = length(replanCobra);
        
        if isempty(replanCobra) 
            break;
        else
            if(firstRun)
                firstRun = false; 
                collCenter = collCenter + sum(distances(1,:)<2);
                numcoll = numcoll + nReplan;
                bad = R2.row;
                % plot(targets(find(distances.dmatrix(:,1)<2),kk),'r.','markersize', 15);
                % plot(targets(bad,kk),'r.','markersize', 15);
            end
            RDSnew = sqrt(rand(nReplan,1) + dA(replanCobra)) .* rRange(replanCobra);
            THTnew = rand(nReplan,1) * 2 * pi;
            targets(replanCobra,kk) =  bench.center(replanCobra) + RDSnew .* exp(1i .* THTnew);
        end
    end
end
R2_ET = toc;
collisions = numcoll;
R2_percentColl = 100* collisions/(numFields * numPos);

numcoll = 0;

%% Rule 3 No collisions on perfect move out. ## Caution this move
%% does not implement our current avoidance strategy. 

currentPosition = bench.center + bench.L1 .* exp(i*bench.tht0) ...
    + bench.L2 .* exp(i*(bench.tht0 + bench.phiIn));

for kk = 1:numFields
    tic
    firstRun = true;
    % this tracks the strategy for each positioner
    usingPrimaryStrategy = true(size(bench.center)); 
    while(true)
        %% calculate the primary trajectory and collisions
        priTrajectories = generateTrajectory(currentPosition,targets(:,kk),...
                                             bench, moveStrategy);
% $$$         priColl = detectCollisions(priTrajectories, bench);
% $$$         priScm = sum(priColl.detected,3); % sum collisions matrix (MxM)
        priColl = detectCollisionsSparse(priTrajectories, bench);
        priScm = accumarray([priColl.row priColl.col], ...
                                  full(sum(priColl.detected,2)), ...
                                  [numPos numPos], [], [], 1);
        priScv = sum(priScm,2);  % sum collisions vector (Mx1)

                           
        %% calculate the alternate trajectory and collisions
        altTrajectories = generateTrajectory(currentPosition,targets(:,kk),...
                                             bench, altmoveStrategy);
% $$$         altColl = detectCollisions(altTrajectories, bench);
% $$$         altScm = sum(altColl.detected,3); % sum collisions matrix
        altColl = detectCollisionsSparse(altTrajectories, bench);
        altScm = accumarray([altColl.row altColl.col], ...
                                  full(sum(altColl.detected,2)), ...
                                  [numPos numPos], [], [], 1);
        altScv = sum(altScm,2); % alt sum collisions vector

        %% for primary trajectories with collisions, replace with
        %% altTrajectory if there is no collision in the
        %% alternative case
        trajectories = priTrajectories;
        scm = priScm;
        scv = priScv;
        repeat_counter = 0;
        last_nCollisionPairs = floor(length(scv)/2);

        while true
            % logical array of cobras that collide
            we_collide = sum(scm,2) | sum(scm,1)';
            switch2alt = we_collide & usingPrimaryStrategy;
            switch2pri = we_collide & ~usingPrimaryStrategy;
            trajectories(switch2alt,:) = altTrajectories(switch2alt,:);
            trajectories(switch2pri,:) = priTrajectories(switch2pri,:);
            
            usingPrimaryStrategy(switch2alt) = false;
            usingPrimaryStrategy(switch2pri) = true;
            
% $$$             fprintf(1,'p->a: %3d\ta->p: %3d\n', length(find(switch2alt)), ...
% $$$                     length(find(switch2pri)));
            
            %% get the collision MxMxN matrix and its projections
% $$$             coll = detectCollisions(trajectories, bench);
% $$$             scm = sum(coll.detected,3); % sum collisions matrix
            coll = detectCollisionsSparse(trajectories, bench);
            scm = accumarray([coll.row coll.col], ...
                             full(sum(coll.detected,2)), ...
                             [numPos numPos], [], [], 1);
            scv = sum(scm,2);  % sum collisions vector
            nCollisionPairs = length(find(we_collide))/2;
            
            if sum(scv) == 0, break; end;
            repeat_counter = repeat_counter + 1;
            fprintf(1,'field %d, replan iteration %d, # collision pairs = %d\n', ...
                    kk, repeat_counter, nCollisionPairs)

            if nCollisionPairs >= last_nCollisionPairs ...
                    & sum(switch2pri) > 0, break ; end;
            last_nCollisionPairs = nCollisionPairs;
        end

        %% plot the trajectory minimum distances.

        minTrajDist = min(coll.minDist,coll.minDist');
        figure(1010)
        cmplx(@scatter,bench.NN.xy,30,nonzeros(minTrajDist),'fill'); 
        hold on;
        plot(bench.NN.xy(find(nonzeros(coll.minDist) < bench.minDist)), ...
             'r.','MarkerSize',30)
        hold off;
        axis equal;
        title('nearest neighbor distances');
        drawnow;

% $$$         figure(1011)
% $$$         histogram(nonzeros(coll.minDist)*2/bench.minDist,1:.25:10);
% $$$         title('distribution of minimum separations over trajectories');

        if sum(scv) == 0
            break;
        else
            fprintf(1,'%d collisions in field %d\n',nnz(scv),kk)
            if(firstRun)
                firstRun = false; 
                numcoll = numcoll + sum(scv > 0); 
                
                if showMoves
                    [rows, cols] = find(scm>0); 
                    for jj = 1:length(rows)
                        rr = rows(jj);
                        cc = cols(jj);
                        if rr < cc
                            figure(jj);
                            showMovement(priTrajectories, bench, priColl, rr, cc);
                            textbp(sprintf('default (primary) trajectories',kk,jj,length(rows)));
                            figure(jj+100);
                            showMovement(altTrajectories, bench, altColl, rr, cc);
                            textbp(sprintf('%s move strategty (Alt strategy)', altmoveStrategy));
                            figure(jj+200);
                            showMovement(trajectories, bench, coll, rr, cc);
                            textbp(sprintf('Final for Field=%d; collision %d/%d',kk,jj,length(rows)));
                        end
                    end
                end
            end 
            
            % decide on which positioner to retarget, assign new target
            disp('Finding new targets...');
            mcc = find(scv == max(scv)); % returns the one who collides with most of the others.
% $$$       figure(9999); plot(targets(mcc(1),kk),'k.','markersize', 18);
            RDSnew = sqrt(rand(1) + dA(mcc(1))) .* rRange(mcc(1));
            THTnew = rand(1) * 2 * pi;
            targets(mcc(1),kk) =  bench.center(mcc(1)) + RDSnew * exp(1i * THTnew);
        end
    end
    if ~mod(kk,100)
        fprintf(1,'Field %d done in %f seconds...\n', kk,toc);
    end;
end

collisions = numcoll/2;
R3_percentColl = 100* collisions/(numFields * numPos);
 
%%
% targetsCenterTP = XY2TP(targets(1,:));
% targetsTP = XY2TP(targets - centers);
% elbowsCenter = targets(1,:) - linkLength * exp(1i* (targetsCenterTP.tht + targetsCenterTP.phi));
% 
% dist = pt2linesegment(targets(2:end,:), repmat(elbowsCenter,6,1), repmat(targets(1,:),6,1));
% collisions = sum(dist<2);
% collisions
% figure(4)
% hold on
% plot(targets(dist(1,:)<2),'r.')


disp('|strategy  | R2 % collisions   | R3 % collisions |')
fprintf(1,'|%-10s|     %4.1f          |       %5.2e  |\n',moveStrategy,...
       R2_percentColl, ...
       R3_percentColl);

clear RDS RDSnew THT THTnew bad centers