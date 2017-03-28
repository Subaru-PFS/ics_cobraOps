function output=simFun(numtrg, cobraLayout, useRealMaps, useRealLinks, varargin)
% collision simulation
% usage: output=simFun(numtrg, cobraLayout, useRealMaps, useRealLinks, varargin)
% varargin: alpha, showMoves
%
% output:
%  targets: Nx1 complex array of target positions
%  Traj   : trajectory structure from realizeTrajectory2
%  Coll   : collision structure from detectCollisionsSparse
%  bench  : bench definition from defineBenchGeometry
% 
% the remaining outputs are meant for code verification.
%
%  minDist: 6Nx1 array of minimum distances between nearest neighbors
%  caats  : "collide at any time" list of nearest neighbors
%  IR1_colliders: not sure
%  IR1          : something to do with interference replan 1

%% Defaults
toggle.info = 'Logical switches';
toggle.showMoves = false;
toggle.showFigures = true;
toggle.SkipTargetReplan = true;
verbosity = 0;
% $$$ thtDIR = 1       ;   % tht & phi have same sense out.
if ~exist('numtrg','var')      , numtrg = 1;          end % number of targets to generate
if ~exist('cobraLayout','var') , cobraLayout = '';    end % use config file bench layout
if ~exist('useRealMaps','var') , useRealMaps  = true; end % use xml data for maps
if ~exist('useRealLinks','var'), useRealLinks = true; end % use xml data for geometry

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
        toggle.showMoves = varargin{jj + 1};
      case 'SkipTargetReplan'
        toggle.SkipTargetReplan = logical(varargin{jj + 1});
      case 'UseThisBench'
        bench = varargin{jj+1};
      case 'verbosity'
        verbosity = varargin{jj+1};
% $$$       case 'thtDIR'
% $$$         thtDIR = sign(varargin{jj + 1});
    end
end

%% performance metrics
PM.info = 'Performance metrics';
PM.total_primary_collisions = 0;
PM.total_collisions = 0;
PM.total_targets    = 0;


if ~exist('bench','var') 
    %% Generate Positioner System
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
    bench = defineBenchGeometry(centers,useRealMaps,useRealLinks);
    clear centers
end
clear useRealLinks useRealMaps cobraLayout varargin;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%     bench defined      %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% reassign bench alpha/beta here (if desired)
if exist('alpha','var')
    bench.alpha = alpha; 
else
    alpha = bench.alpha;
end

numPos = length(bench.center);

if isreal(numtrg)
    TGT_GEN_STRATEGY = 'field'; % uniform over... ('field|'patrol')
else
    TGT_GEN_STRATEGY = 'targetlist';
end

switch TGT_GEN_STRATEGY
  case 'targetlist'
    numFields = 1;
    targets = numtrg(:);
    assignments.tgt = targets; % for partial compatibility with
                               % 'field' case.

  case 'field'
    % us numtrg as density
    numFields = 1;
    targets = complex(zeros(numPos,numFields));
    for kk = 1:numFields
       assignments =  assign_targets(numtrg, bench);
       targets(:,kk) = assignments.tgt;
    end
    PM.R2_percentColl = nan;
  
  case 'patrol'
    %% Test new set of targets using rules:
    numFields = ceil(numtrg / numPos); clear numtrg;
    
    % (R, theta)  coordinate of target
    THT = rand(numPos,numFields)*2*pi;
    RDS = bsxfun(@times, sqrt(bsxfun(@plus, rand(size(THT)), bench.dA)), bench.rRange);
    
    %% positions in patrol area are RDS .* exp(i*THT)
    
    %% Rule 2 No targets closer than 2mm to any line in its target location. 
    
    targets  = bsxfun(@plus, RDS.*exp(1i*THT), bench.center);
    clear RDS THT
    
    numcoll = 0;
    
    for kk = 1:numFields
        toggle.firstRun = true;
        while(true)
            distances = CalcDistanceMatrix(targets(:,kk), bench);
            distances = min(distances, distances');
            [R2.row R2.col] = find(distances > 0 & distances < bench.minDist);
            replanCobra = R2.row(R2.row < R2.col); % only replan half of the cobras.
            nReplan = length(replanCobra);
            
            if isempty(replanCobra) % while loop exit condition
                break;  
            else
                if(toggle.firstRun)
                    toggle.firstRun = false; 
                    numcoll = numcoll + nReplan;
                    bad = R2.row;
                    keyboard;
% $$$                 plot(targets(find(distances.dmatrix(:,1)<2),kk),'r.','markersize', 15);
                    plot(targets(bad,kk),'r.','markersize', 15);
                end
                RDSnew = sqrt(rand(nReplan,1) + bench.dA(replanCobra)) .* bench.rRange(replanCobra);
                THTnew = rand(nReplan,1) * 2 * pi;
                targets(replanCobra,kk) =  bench.center(replanCobra) + RDSnew .* exp(1i .* THTnew);
            end
        end
    end
    PM.R2_ET = toc;
    collisions = numcoll;
    PM.R2_percentColl = 100* collisions/(numFields * numPos);

    clear collisions numcoll replanCobra distances RDSnew THTnew R2
    
  otherwise
    warning(['target strategy is uniform over ''field'' or uniform ' ...
             'over ''patrol''']);
    return;
end
      
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  targets defined, no end-point physical interferences (Rule 2) %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Decide on which hard stop to use.
%  There is a ~25 deg overlap in the tht patrol range.  to prevent having to go after targets just
%  on the wrong side of the HS, set the allowable range from DeltaTheta/2 to MaxThtTravel (~Pi).
%%
% $$$ targetsTP = XY2TP(bsxfun(@minus, targets, bench.center), bench.L1, bench.L2);
% $$$ % delta angles out of same (0) and opposite (1) set ups.
% $$$ dtht0 = mod(targetsTP.tht - bench.tht0, 2*pi);
% $$$ dtht1 = mod(targetsTP.tht - bench.tht1, 2*pi) - 2*pi;
% $$$ tht_overlap = mod(bench.tht1 - bench.tht0,2*pi);
% $$$ max_tht_trvl = pi + tht_overlap;
% $$$ % in deciding which hard stop to use, we allow some overlap both at small angles (1/3 of small
% $$$ % overlap angle) and large angles (pi + overlap)
% $$$ % what follows are logical array cuts on the positioners.
% $$$ ss_allowed =  dtht0 > tht_overlap/3 &  dtht0 < max_tht_trvl;
% $$$ os_allowed = -dtht1 > tht_overlap/3 & -dtht1 < max_tht_trvl; 
% $$$ only_ss = ss_allowed & ~os_allowed;
% $$$ only_os = os_allowed & ~ss_allowed;
% $$$ s_or_os = ss_allowed &  os_allowed;
% $$$ % these take the shorter path, exept when it's too close to zero.
% $$$ use_ss  = s_or_os & (dtht0 < -dtht1) | only_ss;
% $$$ use_os  = s_or_os & (dtht0 > -dtht1) | only_os;

numcoll = 0;

% $$$ currentPosition = bench.home0; % for backwards compatibility.

for kk = 1:numFields
    toggle.firstRun = true;
    % this tracks the strategy for each positioner
    usingPrimaryStrategy = true(size(bench.center)); 
    loopCounter = 0;
    %% Rule 3 No collisions on perfect move out. 
    bench.alpha = 0;
    while(true)
        loopCounter = loopCounter + 1;
        
        % generateTrajectory2 calculates tht and phi vs step.  strategy (decisions on early/late, ss/os)
        % not implemented until realizeTrajectory.
        proto = generateTrajectory2(targets(:,kk), bench);
        % primary strategy is early/late with ss/os picked for fastest route. 
        % realizeTrajectory2 also calculates the late/late os path as the alternate.
        Traj    = realizeTrajectory2(proto,bench);
        Coll    = detectCollisionsSparse(Traj.traj, bench);
        
        repeat_counter = 0;
        last_nCollisions = length(find(Coll.V)); % number of cobras involved in collisions

        if loopCounter == 1 && toggle.firstRun
            PM.total_primary_collisions = PM.total_primary_collisions + last_nCollisions;
        end
        
        last_Traj = Traj;
        last_Coll = Coll;
        
        %% Alter trajectory strategy to minimize collisions on move out
        while true %% exits when "Rule 3" collisions are minimized.
            % index array of cobras that collide
            we_collide = find(last_Coll.V);
            if verbosity, disp('cobra: same-sense?  delta t'); end;
            this_Traj = last_Traj;
            this_usingPrimaryStrategy = usingPrimaryStrategy;
            %            keyboard
            %for jj = 1:length(we_collide)
            for ii = we_collide.';
                if (this_Traj.useP(ii) & this_Traj.ltdiff(ii) > 0) 
                    % true == using same sense tht and going to opp
                    % sense would take more time.
                    fprintf(1,['intractable collision on %4d: use+: ' ...
                    '%d   deltaT: %4d\n'], ii, this_Traj.useP(ii), this_Traj.ltdiff(ii));
                    % in other cases (elseif's here), switch between primary to alternate
                    % strategies depending on what is currently used.
                elseif this_usingPrimaryStrategy(ii)
                    this_Traj.traj(ii,:) = Traj.trajL(ii,:);
                    this_usingPrimaryStrategy(ii) = false;
                elseif ~this_usingPrimaryStrategy(ii)
                    this_Traj.traj(ii,:) = Traj.traj(ii,:);
                    this_usingPrimaryStrategy(ii) = true;
                end
            end

            %% PROBABLY UNNECESSARY TIME SHIFT STUFF
% $$$             [R3.r R3.c] = find(Coll.M); % cobra index of colliders
% $$$             badpair = find(sum(Coll.detected,2))'; % pair index of colliders
% $$$             
% $$$             %% select the cobra from a given pair with the shorter phi move and time shift it
% $$$             phisteps = proto.nphi([R3.r R3.c]);
% $$$             phitime  = proto.lphi([R3.r R3.c]);
% $$$             
% $$$             [tmp indx] = min(phisteps,[],2);
% $$$             % cobras to time shift
% $$$             [indx_short uindx] = unique([R3.r(indx == 1)' R3.c(indx == 2)']);
% $$$             indx_long  = subarray([R3.c(indx == 1)' R3.r(indx == 2)'],uindx);
% $$$ 
% $$$             phiDT = sparse(1,2394);
% $$$             phiDT(indx_long) = proto.lphi(indx_short);
% $$$             keyboard
% $$$             T = modifyTrajectory(Traj,proto,bench,phiDT);

            this_Coll = detectCollisionsSparse(this_Traj.traj,bench);
            still_collide = find(this_Coll.V);
            nCollisions = length(still_collide);
            repeat_counter = repeat_counter + 1;
                                     
            if (verbosity>0) 
                fprintf(1,'Field %d InterferenceReplan %02d #StrategyReplan%d  #Coll = %d -> %2d\n',...
                        kk, loopCounter, repeat_counter, length(we_collide), ...
                        nCollisions);
            end

            %%% temporary::: show all collisions
            [rows, cols] = find(this_Coll.M);
            if toggle.showMoves & false
                for jj = 1:length(rows)
                    rr = rows(jj);
                    cc = cols(jj);
                    if rr < cc
                        figure(1000+rr);
                        showMovementNN(Traj.traj, bench, Coll, rr, targets);
                        textbp(sprintf('default (primary) trajectories',kk,jj,length(rows)));
                        figure(10000+rr);
                        showMovementNN(this_Traj.traj, bench, this_Coll, rr, targets);
                        textbp(sprintf('alt? move strategy'));
                    end
                end
            end
            %%%%

            if nCollisions >= last_nCollisions % & sum(switch2pri) > 0
                % in here means that switching stragetgies did not improve the situation
                break ; 
            end
            usingPrimaryStrategy = this_usingPrimaryStrategy;
            last_nCollisions = nCollisions;
            last_Traj = this_Traj;
            last_Coll = this_Coll;
            if nCollisions == 0, break; end;
        end % R3, inner while loop -solve collisions via alternate trajectories.
        %% use last_Traj and last_Coll as the working trajectory.
        Traj = last_Traj;
        Coll = last_Coll;
        nCollisions = last_nCollisions;

        clear last_Traj last_Coll this_Traj this_Traj last_nCollisions
        
        minTrajDist = min(Coll.minDist,Coll.minDist');
        %% plot the trajectory minimum distances.
% $$$         if true
% $$$             figure(1010)
% $$$             try
% $$$                 cmplx(@scatter,bench.NN.xy,30,nonzeros(minTrajDist),'fill'); 
% $$$             catch
% $$$                 indx = sub2ind(size(minTrajDist), bench.NN.row, bench.NN.col);
% $$$                 zind = indx(find(minTrajDist(indx) == 0))
% $$$                 [r c] = ind2sub(size(minTrajDist),zind)
% $$$                 cmplx(@scatter,bench.NN.xy,30,...
% $$$                       minTrajDist(sub2ind(size(minTrajDist),bench.NN.row,bench.NN.col)),'fill'); 
% $$$                 disp('look at minTrajDist');
% $$$                 % this problem originates in generateTrajectory.
% $$$                 keyboard;
% $$$             end;
% $$$             hold on;
% $$$             plot(bench.NN.xy(find(nonzeros(Coll.minDist) < bench.minDist)), ...
% $$$                  'r.','MarkerSize',30)
% $$$             hold off;
% $$$             axis equal;
% $$$             title('nearest neighbor distances');
% $$$             colorbar; 
% $$$             colormap([1 0 0;...
% $$$                       1 0 0;...
% $$$                       1 1 0;...
% $$$                       1 1 0;...
% $$$                       0 0 1;...
% $$$                       0 0 1;...
% $$$                       0 0 1;...
% $$$                       0 0 1;...
% $$$                       0 0 1;...
% $$$                       0 0 1;...
% $$$                      ]);
% $$$             caxis([0,10]);
% $$$             drawnow;
% $$$         end 
        %% end of plot min trajectory distances

        minDist(:,loopCounter) = nonzeros(minTrajDist);

        if toggle.showFigures & false
            figure(1011)
            %% asymmetric distances
            % histogram(nonzeros(Coll.minDist)*2/bench.minDist,0:.25:10);
            %% symmetric distances
            histogram(minDist*2/bench.minDist,0:.25:10);
            xlabel('minimal fiber to arm distance [mm]');
            title('distribution of minimum separations over trajectories');
        end

        %% collect/summarize failures
        if toggle.firstRun
            PM.total_collisions = PM.total_collisions + nCollisions;
            PM.total_targets    = PM.total_targets + numPos;
            IR1_colliders    = find(we_collide);
            IR1_tp           = XY2TP(targets - bench.center, bench.L1, bench.L2);
            [r c] = find(bench.nnMap(IR1_colliders, IR1_colliders));
            keep = r<c;
            r = IR1_colliders(r(keep));
            c = IR1_colliders(c(keep));
            IR1.ax  = angle(bench.center(c) - bench.center(r))';
            IR1.tht = [IR1_tp.tht(r) IR1_tp.tht(c)]';
            IR1.phi = [IR1_tp.phi(r) IR1_tp.phi(c)]';
            IR1.coll = Coll;
            IR1.traj = Traj.traj;
            IR1.UPS  = usingPrimaryStrategy;
            clear r c keep
        end

        if nCollisions == 0 || toggle.SkipTargetReplan
            break;
        else % there are still collisions
            %%% FINDING NEW TARGETS FOR COBRAS
            [rows, cols] = find(Coll.M); 
            if (toggle.firstRun)
                toggle.firstRun = false; 
                numcoll = numcoll + nCollisions; 

                %% [2016 07 18:  why do I need to plot these?
% $$$                 if toggle.showMoves
% $$$                     for jj = 1:length(rows)
% $$$                         rr = rows(jj);
% $$$                         cc = cols(jj);
% $$$                         if rr < cc
% $$$                             figure(12);
% $$$                             showMovementNN(  primary.traj, bench, priColl, rr, targets(:,kk));
% $$$                             try
% $$$                                 textbp(sprintf('default (primary) trajectories',kk,jj,length(rows)));
% $$$                             catch
% $$$                             end
% $$$                             figure(12+100);
% $$$                             showMovementNN(alternate.traj, bench, altColl, rr, targets(:,kk));
% $$$                             try
% $$$                                 textbp(sprintf('%s move strategy (Alt strategy)', altmoveStrategy));
% $$$                             catch
% $$$                             end
% $$$                             figure(12+200);
% $$$                             showMovementNN(  trajectories, bench,    coll, rr, targets(:,kk)); 
% $$$                             try
% $$$                                 textbp(sprintf('Final for Field=%d; collision %d/%d',kk,jj,length(rows)));
% $$$                             catch
% $$$                             end
% $$$                         end
% $$$                     end
% $$$                 end
            end 
            
            % decide on which positioner to retarget, assign new target

            replanCobras = [];
            for jj = 1:length(rows)
                rr = rows(jj);
                cc = cols(jj);
                disp('Need to fix things here to proceed with target replanning (simFun L449)')
                keyboard
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
                RDSnew = sqrt(rand(size(replanCobras)) + bench.dA(replanCobras)) .* bench.rRange(replanCobras);
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
            end % target replan loop
        end % target replan if branch
    end % outer while loop: bad targets replaced.  continuing.

    %% trajetories that are NOT using the primary strategy are, by
    %% definition, using negative theta (opp sense) and late theta
    Traj.useP(~usingPrimaryStrategy) = false;
    Traj.useL(~usingPrimaryStrategy) = true;

    if ~mod(kk,100)
        fprintf(1,'Field %d done in %f seconds...\n', kk,toc);
    end;
    bench.alpha = alpha;
end % loop over fields.

PM.R3_percentColl = 100* nCollisions/(numFields * numPos);

if (verbosity > 0)
    disp('|strategy  | R2 % collisions   | R3 % collisions |')
    fprintf(1,'|%-10s|     %4.1f          |       %5.2e  |\n','dual HS',...
            PM.R2_percentColl, ...
            PM.R3_percentColl);
end

caats = find(min(minDist,[],2) < bench.minDist); % index of fibers that have Collided At Any Time
clate = find(min(minDist(:,5:end),[],2) < bench.minDist); % index of fibers that have Collided LATE
cstfr = caats(find(minDist(caats,1) > 4)); % indx of fibs that Collide any time but STart FaR

%%
% targetsCenterTP = XY2TP(targets(1,:));
TP = XY2TP(targets - bench.center, bench.L1, bench.L2);
% elbowsCenter = targets(1,:) - linkLength * exp(1i* (targetsCenterTP.tht + targetsCenterTP.phi));
% 
% dist = pt2linesegment(targets(2:end,:), repmat(elbowsCenter,6,1), repmat(targets(1,:),6,1));
% collisions = sum(dist<2);
% collisions
% figure(4)
% hold on
% plot(targets(dist(1,:)<2),'r.')

if toggle.showFigures
% $$$     figure(380);
% $$$     %% patrol areas and targets
% $$$     plot(assignments.tgt,'k.'); hold on;
% $$$     plot(assignments.rem,'k.');
% $$$     cmplx(@plotcircle, bench.center, bench.rMax, 'b-');
% $$$     hold off;
% $$$     axis equal;

    figure(381);
    %% targets and collisions
    cobra_arms = [bench.center ...
                  bench.center + bench.L1.*exp(i*TP.tht) ...
                  targets].'; % cobra [center; elbow; fiber] array (3 x ncobras)
    all_coll_mindist = nonzeros(Coll.minDist);
    colliders  = find(all_coll_mindist <= bench.minDist);
    nearcolliders = find(all_coll_mindist < bench.minDist*1.5 & ...
                         all_coll_mindist > bench.minDist); %% watch out -- HARD CODE
    plot(cobra_arms,'ko-');
    hold on;
    if strcmp(TGT_GEN_STRATEGY, 'field')
        plot(assignments.rem,'bx'); % unassigned targets from full list
        plot(targets(assignments.isassigned),'bo','markerface','g');
        plot(targets(~assignments.isassigned),'bo','markerface','k');
    else
        plot(targets,'bo','markerface','g');
    end;
    plot(bench.NN.xy(nearcolliders),'yo','MarkerFace','y');
    plot(bench.NN.xy(colliders),'ro','MarkerFace','r');
    plot(Traj.traj.','k--');
    hold off;
    cmplx(@plotcircle, bench.center, bench.rMax, 'k:');
    axis equal;
    %    set(gca,'Xdir','reverse','Ydir','reverse');
    if length(colliders) > 0
        cmplx(@text, bench.center(bench.NN.row(colliders)), cellstr(num2str(bench.NN.row(colliders), '%05d')));
        
        %% show cobras that collide -- VERY SLOW
        if toggle.showMoves
            for pid = unique(min([bench.NN.col(colliders) bench.NN.row(colliders)],[],2))'
                hold on;
                showMovementNN(Traj.traj, bench, Coll, pid, targets);
                hold off;
            end
        end
    end
end

clear cobra_arms all_coll_mindist colliders nearcolliders

% $$$ figure(382);
% $$$ plot(bench.NN.xy,'g.'); axis equal; hold on;
% $$$ plot(bench.NN.xy(caats),'b.','MarkerSize',30);
% $$$ %plot(bench.NN.xy(clate),'r.','MarkerSize',20);
% $$$ %plot(bench.NN.xy(cstfr),'kx','MarkerSize',20);
% $$$ hold off;
% $$$ try
% $$$     legend('inter-Cobra locations','collides')%,'collides late','strt w/ big clearance');
% $$$     cmplx(@text,bench.NN.xy(caats),cellstr(num2str(caats)));
% $$$ catch
% $$$ end
% $$$ title('Map of all collisions over replan history');

if ~toggle.SkipTargetReplan
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
end

output = packstruct(...
    targets ... % desired target locations
    ,Traj ... % final trajectories
    ,Coll ... % collision structure
    ,bench ... % bench geometry structure
    ,minDist ... % minimum distance between adjacent cobras on each loop
    ,caats ... % standard analysis outputs for PHM
    ,assignments...
    );
% $$$     ,nSteps ... % number of steps requested for each motor
% $$$     ,IR1_colliders... % cobras that collide at end of first interference replan
% $$$     ,IR1...

% $$$ fprintf(1,'%f\n',PM.total_collisions/PM.total_targets)
