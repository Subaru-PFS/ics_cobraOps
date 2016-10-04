function output=simFun(numtrg, cobraLayout, useRealMaps, useRealLinks, varargin)
% This version of the collision simulation produces colliding sets. It
% needs a good balance between maxNumColl (maximum number of collisions)
% and the attempts to generate these. Count your cobras. 
maxNumColl = 18;
attempts = 60;
% usage: output=simFun(numtrg, cobraLayout, useRealMaps, useRealLinks, varargin)
% varargin: alpha, showMoves, thtDIR

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
    end
end

%% Defaults
if ~exist('numtrg','var')      , numtrg = 1000;       end % number of targets to generate
if ~exist('showMoves','var')   , showMoves = false;   end % Show single collision movements
if ~exist('thtDIR','var')      , thtDIR = 1       ;   end % tht & phi have same sense out.
if ~exist('useRealMaps','var') , useRealMaps  = true; end % use xml data for maps
if ~exist('useRealLinks','var'), useRealLinks = true; end % use xml data for geometry

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
minDist = [];
%% Rule 3 No collisions on perfect move out.

currentPosition = bench.center + bench.L1 .* exp(i*bench.tht0) ...
    + bench.L2 .* exp(i*(bench.tht0 + bench.phiIn));

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
        
%         fprintf(1,'Field %d Gen Traj %02d #  ',...
%             kk, loopCounter);
        
        
        %% for primary trajectories with collisions, replace with
        %% altTrajectory if there is no collision in the
        %% alternative case
        trajectories = primary.traj;
        scm = priScm;
        scv = priScv;
        
        %% plot the trajectory minimum distances.
        
        if sum(scv > 0) > maxNumColl  
            
            numColl = sum(scv > 0);
            fprintf(1,'Done, reached: %d collisions', numColl + 1 - 1);
            
            break;
        elseif (loopCounter > attempts)
            numColl = sum(scv > 0);
            fprintf(1,'Tried so hard, couldnt reach more than: %d collisions', numColl + 1 - 1);
            break;
        else
            [rows, cols] = find(scm);
            if(firstRun)
                firstRun = false;
                numcoll = numcoll + sum(scv > 0);
            end
            
            replanCobras = [];
            repeat_counter = 0;
            remTarg = ones(size(scm));
            remTarg(cols,:) = 0;
            remTarg(:,rows) = 0;
            %
            [rrow, rcol] = find(remTarg > 0);
            replanCobras = unique(rrow);
            while size(replanCobras, 1) > 2;
                repeat_counter = repeat_counter + 1;
%                 fprintf(1,'Field %d #Target__Replan%2d  #NotColliding = %2d\n',...
%                     kk, loopCounter, length(replanCobras));
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

caats = find(min(minDist,[],2) < bench.minDist); % index of fibers that have Collided At Any Time
delta = bsxfun(@minus,[targets(:,end) - trajectories(:,end)], bench.center); % position errors in local cobra coordinates.



%% reconstruct nSteps structure with the actual trajectory
%% information
nSteps = rmfield(primary,'traj');
%nSteps.dtht(usingPrimaryStrategy) = primary.dtht(usingPrimaryStrategy);
%nSteps.dphi(usingPrimaryStrategy) = primary.dphi(usingPrimaryStrategy);


output = packstruct(...
    targets ... % desired target locations
    ,trajectories ... % final trajectories
    ,bench ... % bench geometry structure
    ,nSteps...
    ,minDist ... % minimum distance between adjacent cobras on each loop
    ,caats,delta ... % standard analysis outputs for PHM
    );