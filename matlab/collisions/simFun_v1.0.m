function output=simFun(numtrg, cobraLayout, useRealMaps, useRealLinks, varargin)
% collision simulation
% usage: output=simFun(numtrg, cobraLayout, useRealMaps, useRealLinks, varargin)
% varargin: alpha, showMoves, thtDIR
%
% this is v 1.0.  single hard stop implemented.  some code for 2
% hard stop setup is here, but not propagated into Rule #3.
%


%% parse varargin
if mod(length(varargin),2) ~= 0
    disp('Varagin should be a list of alternating parameter names and values');
    disp('Exiting...');
    return;
end
for jj=1:2:length(varargin)
    %% need to put some input parsing code here.
    switch varargin{jj}
      case 'alpha'
        alpha = varargin{jj + 1};
      case 'showMoves'
        showMoves = varargin{jj + 1};
      case 'thtDIR'
        thtDIR = sign(varargin{jj + 1});
      case 'SkipTargetReplan'
        SkipTargetReplan = logical(varargin{jj + 1});
    end
end

%% Defaults
if ~exist('numtrg','var')      , numtrg = 1000;       end % number of targets to generate
if ~exist('showMoves','var')   , showMoves = false;   end % Show single collision movements
if ~exist('thtDIR','var')      , thtDIR = 1       ;   end % tht & phi have same sense out.
if ~exist('useRealMaps','var') , useRealMaps  = true; end % use xml data for maps
if ~exist('useRealLinks','var'), useRealLinks = true; end % use xml data for geometry
if ~exist('SkipTargetReplan','var')
    SkipTargetReplan = false;
end;

%% performance metrics
total_primary_collisions = 0;
total_collisions = 0;
total_targets    = 0;


%% Specify move strategy
%moveStrategy = 'earlyEarly';
%altmoveStrategy = 'earlyLate';
moveStrategy = 'earlyLate'; % aka late, or L hook
altmoveStrategy = 'lateLate'; % aka radial, or R hook

%% Generate Positioner System
%% Q: why is the 7 pos system worse than the 14 rail system?

switch cobraLayout
  case {'','none'}
    % use config geometry
    centers = [];
    useRealMaps = true;
  case 'hex'
    % 7 positioners in hex pattern      (R3 collisions = 3.0e-4 w/ 100k targets)
    centers = [0 8 * exp(1i * (0:5) * pi / 3)].';
  case 'line'
    % Line of cobras                    (R3 collisions = 4.5e-5 (fraction) w/ 1M targets)
    centers = complex([0:8:(26*8)],0)' + i;
  case 'rails'
    % N hard coded rails                (R3 collisions = 9.9e-5 w/ 16k targets)
    centers = getCentersRails(14);
  case 'full'
    % full PFI bench                    (R3 collisions = 1.7e-4 (fraction) w/ 1M targets)
    centers = getCentersRails(14);
    centers = bsxfun(@times, [1 exp(i*2*pi/3) exp(i*4*pi/3)], centers);
    centers = reshape(centers,[],1);
  otherwise
    disp('Valid cobra layouts: ''none'',''line'',''rails'',''full''');
    disp('Exiting...');
    return;
end
    
%% Define a structure with the bench geometry

bench = defineBenchGeometry(centers,useRealMaps,useRealLinks,thtDIR);

% reassign bench alpha/beta here (if desired)
if exist('alpha','var'), bench.alpha=alpha; end;

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


%% Decide on which hard stop to use.
%  There is a ~25 deg overlap in the tht patrol range.  to prevent having to go after targets just
%  on the wrong side of the HS, set the allowable range from DeltaTheta/2 to MaxThtTravel (~Pi).
%%
targetsTP = XY2TP(bsxfun(@minus, targets, bench.center), bench.L1, bench.L2);
% delta angles out of same (0) and opposite (1) set ups.
dtht0 = mod(targetsTP.tht - bench.tht0, 2*pi);
dtht1 = mod(targetsTP.tht - bench.tht1, 2*pi) - 2*pi;
tht_overlap = mod(bench.tht1 - bench.tht0,2*pi);
max_tht_trvl = pi + tht_overlap;
% in deciding which hard stop to use, we allow some overlap both at small angles (1/3 of small
% overlap angle) and large angles (pi + overlap)
% what follows are logical array cuts on the positioners.
ss_allowed =  dtht0 > tht_overlap/3 &  dtht0 < max_tht_trvl;
os_allowed = -dtht1 > tht_overlap/3 & -dtht1 < max_tht_trvl; 
only_ss = ss_allowed & ~os_allowed;
only_os = os_allowed & ~ss_allowed;
s_or_os = ss_allowed &  os_allowed;
% these take the shorter path, exept when it's too close to zero.
use_ss  = s_or_os & (dtht0 < -dtht1) | only_ss;
use_os  = s_or_os & (dtht0 > -dtht1) | only_os;

numcoll = 0;

%% Rule 3 No collisions on perfect move out. 

ss_home = bench.center + bench.L1 .* exp(i*bench.tht0) + bench.L2 .* exp(i*(bench.tht0 + bench.phiIn));
os_home = bench.center + bench.L1 .* exp(i*bench.tht1) + bench.L2 .* exp(i*(bench.tht1 + bench.phiIn));
currentPosition = ss_home; % for backwards compatibility.

keyboard;

for kk = 1:numFields
    tic
    firstRun = true;
    % this tracks the strategy for each positioner
    usingPrimaryStrategy = true(size(bench.center)); 
    loopCounter = 0;
    while(true)
        loopCounter = loopCounter + 1;
        %% calculate the primary trajectory and collisions
        primary = generateTrajectory(currentPosition,targets(:,kk),...
                                     bench, moveStrategy);
        priColl = detectCollisionsSparse(primary.traj, bench);
        priScm = accumarray([priColl.row priColl.col], ...
                                  full(sum(priColl.detected,2)), ...
                                  [numPos numPos], [], [], 1);
        priScv = sum(priScm,2);  % sum collisions vector (Mx1)

                           
        %% calculate the alternate trajectory and collisions
        alternate = generateTrajectory(currentPosition,targets(:,kk),...
                                             bench, altmoveStrategy);
        altColl = detectCollisionsSparse(alternate.traj, bench);
        altScm = accumarray([altColl.row altColl.col], ...
                                  full(sum(altColl.detected,2)), ...
                                  [numPos numPos], [], [], 1);
        altScv = sum(altScm,2); % alt sum collisions vector

        %% for primary trajectories with collisions, replace with
        %% altTrajectory if there is no collision in the
        %% alternative case
        trajectories = primary.traj;
        scm = priScm;
        scv = priScv;
        repeat_counter = 0;
        last_nCollisionPairs = floor(length(find(scv))/2);

        if loopCounter == 1 && firstRun
            total_primary_collisions = total_primary_collisions + last_nCollisionPairs;
        end

        while true
            % logical array of cobras that collide
            we_collide = sum(scm,2) | sum(scm,1)';
            switch2alt = we_collide & usingPrimaryStrategy;
            switch2pri = we_collide & ~usingPrimaryStrategy;
            trajectories(switch2alt,:) = alternate.traj(switch2alt,:);
            trajectories(switch2pri,:) = primary.traj(switch2pri,:);
            
            usingPrimaryStrategy(switch2alt) = false;
            usingPrimaryStrategy(switch2pri) = true;
            
% $$$             fprintf(1,'p->a: %3d\ta->p: %3d\n', length(find(switch2alt)), ...
% $$$                     length(find(switch2pri)));
            
            %% get the collision MxMxN matrix and its projections
            coll = detectCollisionsSparse(trajectories, bench);
            scm = accumarray([coll.row coll.col], ...
                             full(sum(coll.detected,2)), ...
                             [numPos numPos], [], [], 1);
            scv = sum(scm,2);  % sum collisions vector
            nCollisionPairs = length(find(we_collide))/2;
            
            if sum(scv) == 0, break; end;
            repeat_counter = repeat_counter + 1;
% $$$             fprintf(1,'field %d, strategy replan iteration %d, # collision pairs = %d\n', ...
% $$$                     kk, repeat_counter, nCollisionPairs)
                                     
            fprintf(1,'Field %d InterferenceReplan %02d #StrategyReplan%d  #Coll = %2d\n',...
                    kk, loopCounter, repeat_counter, ceil(nCollisionPairs))

            %%% temporary::: show all collisions
            [rows, cols] = find(scm);
            if showMoves
                for jj = 1:length(rows)
                    rr = rows(jj);
                    cc = cols(jj);
                    if rr < cc
                        figure(jj);
                        showMovementNN(primary.traj, bench, priColl, rr, targets);
                        textbp(sprintf('default (primary) trajectories',kk,jj,length(rows)));
                        figure(jj+100);
                        showMovementNN(alternate.traj, bench, altColl, rr, targets);
                        textbp(sprintf('%s move strategy (Alt strategy)', altmoveStrategy));
                        figure(jj+200);
                        showMovementNN(   trajectories, bench,    coll, rr, targets); 
                        textbp(sprintf('Final for Field=%d; collision %d/%d',kk,jj,length(rows)));
                        keyboard;
                    end
                end
            end
            %%%%


            if nCollisionPairs >= last_nCollisionPairs ...
                    & sum(switch2pri) > 0
                % in here means that switching stragetgies made
                % things worse.
                break ; 
            end
            last_nCollisionPairs = nCollisionPairs;
        end % R3, inner while loop -solve collisions via alternate trajectories.

        
        %% plot the trajectory minimum distances.
        if false  
            minTrajDist = min(coll.minDist,coll.minDist');
            figure(1010)
            try
                cmplx(@scatter,bench.NN.xy,30,nonzeros(minTrajDist),'fill'); 
            catch
                indx = sub2ind(size(minTrajDist), bench.NN.row, bench.NN.col);
                zind = indx(find(minTrajDist(indx) == 0))
                [r c] = ind2sub(size(minTrajDist),zind)
                cmplx(@scatter,bench.NN.xy,30,...
                      minTrajDist(sub2ind(size(minTrajDist),bench.NN.row,bench.NN.col)),'fill'); 
                disp('look at minTrajDist');
                % this problem originates in generateTrajectory.
                keyboard;
            end;
            hold on;
            plot(bench.NN.xy(find(nonzeros(coll.minDist) < bench.minDist)), ...
                 'r.','MarkerSize',30)
            hold off;
            axis equal;
            title('nearest neighbor distances');
            colorbar; 
            colormap([1 0 0;...
                      1 0 0;...
                      1 1 0;...
                      1 1 0;...
                      0 0 1;...
                      0 0 1;...
                      0 0 1;...
                      0 0 1;...
                      0 0 1;...
                      0 0 1;...
                     ]);
            caxis([0,10]);
            drawnow;
        end 
        %% end of plot min trajectory distances

% $$$         fprintf(1,'# steps requested = %d\n',ceil(nSteps.max));

% $$$         figure(1011)
% $$$         histogram(nonzeros(coll.minDist)*2/bench.minDist,1:.25:10);
% $$$         title('distribution of minimum separations over trajectories');

        minDist(:,loopCounter) = nonzeros(coll.minDist);
        
        %% collect/summarize failures
        if firstRun
            total_collisions = total_collisions + last_nCollisionPairs;
            total_targets    = total_targets + length(centers);
            IR1_colliders    = find(we_collide);
            IR1_tp           = XY2TP(targets - bench.center, bench.L1, bench.L2);
            [r c] = find(bench.nnMap(IR1_colliders, IR1_colliders));
            keep = r<c;
            r = IR1_colliders(r(keep));
            c = IR1_colliders(c(keep));
            IR1.ax  = angle(bench.center(c) - bench.center(r))';
            IR1.tht = [IR1_tp.tht(r) IR1_tp.tht(c)]';
            IR1.phi = [IR1_tp.phi(r) IR1_tp.phi(c)]';
            IR1.coll = coll;
            IR1.traj = trajectories;
            IR1.UPS  = usingPrimaryStrategy;
            clear r c keep
        end
        if SkipTargetReplan
            break 
        end
        if sum(scv) == 0 % no collisions
            break;
        else % there are collisions still
% $$$             fprintf(1,'%d collisions in field %d\n',nnz(scv),kk)
            [rows, cols] = find(scm); 
            if(firstRun)
                firstRun = false; 
                numcoll = numcoll + sum(scv > 0); 
                
                if showMoves
                    for jj = 1:length(rows)
                        rr = rows(jj);
                        cc = cols(jj);
                        if rr < cc
                            figure(12);
                            showMovementNN(  primary.traj, bench, priColl, rr, targets(:,kk));
                            try
                                textbp(sprintf('default (primary) trajectories',kk,jj,length(rows)));
                            catch
                            end
                            figure(12+100);
                            showMovementNN(alternate.traj, bench, altColl, rr, targets(:,kk));
                            try
                                textbp(sprintf('%s move strategy (Alt strategy)', altmoveStrategy));
                            catch
                            end
                            figure(12+200);
                            showMovementNN(  trajectories, bench,    coll, rr, targets(:,kk)); 
                            try
                                textbp(sprintf('Final for Field=%d; collision %d/%d',kk,jj,length(rows)));
                            catch
                            end
                        end
                    end
                end
            end 
            
            % decide on which positioner to retarget, assign new target

            replanCobras = [];
            for jj = 1:length(rows)
                rr = rows(jj);
                cc = cols(jj);
                % from each pair, take the target with the longest theta moves
                if primary.ntht(rr) > primary.ntht(cc)
                    replanCobras = [replanCobras; rr];
                else
                    replanCobras = [replanCobras; cc];
                end
                replanCobras = unique(replanCobras);
            end
            repeat_counter = 0;
            while ~isempty(replanCobras)
                repeat_counter = repeat_counter + 1;
                fprintf(1,'Field %d InterferenceReplan %02d #Target__Replan%d  #Coll = %2d\n',...
                        kk, loopCounter, repeat_counter, length(replanCobras))
                RDSnew = sqrt(rand(size(replanCobras)) + dA(replanCobras)) .* rRange(replanCobras);
                THTnew = rand(size(replanCobras)) * 2 * pi;
                targets(replanCobras,kk) =  bench.center(replanCobras) + RDSnew .* exp(1i * THTnew);
                distances = CalcDistanceMatrix(targets(:,kk), bench);
                distances = min(distances,distances.');
                retry = false(size(replanCobras));
                for indx = 1:length(replanCobras)
                    if min(nonzeros(distances(replanCobras(indx),:))) < bench.minDist
                        retry(indx) = true;
                    end
                end
                replanCobras = replanCobras(retry);
            end
        end
    end
    if ~mod(kk,100)
        fprintf(1,'Field %d done in %f seconds...\n', kk,toc);
    end;
end % loop over fields.

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

caats = find(min(minDist,[],2) < bench.minDist); % index of fibers that have Collided At Any Time
clate = find(min(minDist(:,5:end),[],2) < bench.minDist); % index of fibers that have Collided LATE
cstfr = caats(find(minDist(caats,1) > 4)); % indx of fibs that Collide any time but STart FaR
delta = bsxfun(@minus,[targets(:,end) - trajectories(:,end)], bench.center); % position errors in local cobra coordinates.

figure(382);
plot(bench.NN.xy,'g.'); axis equal; hold on;
plot(bench.NN.xy(caats),'b.','MarkerSize',30);
%plot(bench.NN.xy(clate),'r.','MarkerSize',20);
%plot(bench.NN.xy(cstfr),'kx','MarkerSize',20);
hold off;
legend('inter-Cobra locations','collides')%,'collides late','strt w/ big clearance');
try
cmplx(@text,bench.NN.xy(caats),cellstr(num2str(caats)));
catch
end
title('Map of all collisions over replan history');


% for colliders, plot minDist history sorted by first replan iteration with collision
[rr cc] = find(minDist < bench.minDist);
[ur ir] = unique(rr);  % unique collisions
[sc ic] = sort(cc(ir));% sort by first-collision order
isu = ir(ic);          % index of sort(unique(collisions))

figure(383);
imagesc(minDist(rr(isu),:));
title(['minimum distance histories for colliding pairs, sorted by ' ...
       'replan index']);
xlabel('target replan index');
ylabel('positioner [arb indx]');
colorbar;

[sd id] = sort(minDist(rr(ir),1));
idsu = ir(id);

figure(384);
imagesc(minDist(rr(idsu),:));
title(['minimum distance histories for colliding pairs, sorted by ' ...
       'initial nearest distance']);
xlabel('target replan index');
ylabel('positioner [arb indx]');
colorbar;

%% reconstruct nSteps structure with the actual trajectory
%% information
nSteps = rmfield(alternate,'traj');
nSteps.dtht(usingPrimaryStrategy) = primary.dtht(usingPrimaryStrategy);
nSteps.dphi(usingPrimaryStrategy) = primary.dphi(usingPrimaryStrategy);

output = packstruct(...
    targets ... % desired target locations
    ,trajectories ... % final trajectories
    ,coll ... % collision structure
    ,bench ... % bench geometry structure
    ,nSteps ... % number of steps requested for each motor
    ,minDist ... % minimum distance between adjacent cobras on each loop
    ,caats,delta ... % standard analysis outputs for PHM
    ,IR1_colliders... % cobras that collide at end of first interference replan
    ,IR1...
    );

fprintf(1,'%f\n',total_collisions/total_targets)