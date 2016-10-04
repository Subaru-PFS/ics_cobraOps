function output = detectCollisionsSparse(Trajectory, geom)
% INPUTS: bench geometry and NN map (sparse MxM matrix)
%         trajectories              (MxN array)
% OUTPUT: return a logical array indicating positioners that violate
%         the minimum fiber to neighboring phi-arm separation

  nCobra = size(Trajectory,1); %[M]
  nSteps = size(Trajectory,2); %[N]
  
  [row col] = find(geom.nnMap); % [6M]
  rcvector  = find(geom.nnMap); % [6M] rcvector for first position in
                                % trajectory
  %% need to calculate matrix indexes for all relevant calcs
  rctvector = reshape( bsxfun( @plus,...
                              rcvector,...
                              (0:nSteps-1) * nCobra^2 ),...
                       [],1);    % [6MN] in row, col, traj ordering
  trcvector = reshape( bsxfun( @plus,...
                               (rcvector-1) * nSteps,...
                               1:nSteps), ...
                       [],1);   % [6MN] in traj, row, col ordering
  % compute the elbow position for each target
  TrajectoryTP = XY2TP(bsxfun(@minus,Trajectory,geom.center), geom.L1, geom.L2);
  Elbows = Trajectory - bsxfun(@times, ...
                               exp(1i* (TrajectoryTP.tht + TrajectoryTP.phi)),...
                               geom.L2); % both are [M x N]

  % fiber and elbow trajectories
  FIB1 = Trajectory(row,:); % [6M x N]
  FIB2 = Trajectory(col,:); % [6M x N]
  ELB2 = Elbows(col,:);     % [6M x N]
  
  % calculate FIB1 - ARM2 distances
  [distance distanceType] = pt2linesegment(FIB1, ELB2, FIB2); % [6M x N]

% $$$   tic
% $$$   detected = (distance < geom.minDist);
% $$$   type     = distanceType .* detected;
% $$$   toc
  
  detected = sparse(distance < geom.minDist);
  type     = distanceType .* detected;
  % type = 1: fiber hits elbow
  % type = 2: fiber hits arm (between elbow and fib)
  % type = 3: fiber hits fiber


  minDist = accumarray([row col], min(distance,[],2), [nCobra nCobra],[],[],1);
  
  %% johannes's sum of collisions matrix and vector
  %  time bins with collisions for each pair (M) or each positioner (V)
  M = accumarray([row col], full(sum(detected,2)), [nCobra nCobra],[],[],1);
  V = sum(M,2) + sum(M,1)';

% $$$   indx_detected = find(sum(detected,2));
% $$$   rc2indx = sparse(row(indx_detected),col(indx_detected),indx_detected);
  rcindx = sparse(row,col,(1:length(row))');
  
  output = packstruct(row, col, rcindx, detected, type, minDist,M,V);
  
  %% verification plots
  if false

    for jj = 1:nSteps
      clf;
      subplot(121)
      cmplx(@plotcircle,geom.center,geom.L1 + geom.L2,'r');
      hold on;
      plot(Trajectory(:,jj),'o');
      plot([Elbows(:,jj) Trajectory(:,jj)].');
      hold off;
      
      dmap = distanceTENSOR(:,:,jj);
      dmap(distanceTENSOR(:,:,jj) > 2) = 1;
      dmap(distanceTENSOR(:,:,jj) < 2) = 2;
      dmap(isnan(dmap)) = 0;
      
      subplot(122)
      imagesc(dmap);
      colorbar;
      caxis([0,2])
      axis equal
      title(sprintf('%d',jj));
      pause(1);
    end
  end
    
  