function output = detectCollisions(Trajectory, geom)
% given Centers (Mx1 vector) and a set of fiber position trajectories
% (MxN array), return a logical array indicating positioners that
% violate the minimum fiber to neighboring phi-arm separation

  nCobra = size(Trajectory,1);
  nSteps = size(Trajectory,2);
  
  [row col] = find(geom.nnMap);
  rcvector  = find(geom.nnMap); % rcvector for first position in
                                % trajectory
  %% need to calculate matrix indexes for all relevant calcs
  rcvector = reshape( bsxfun( @plus,...
                              rcvector,...
                              (0:nSteps-1) * nCobra^2 ),...
                      [],1);
  
  % compute the elbow position for each target
  TrajectoryTP = XY2TP(bsxfun(@minus,Trajectory,geom.center), geom.L1, geom.L2);
  Elbows = Trajectory - bsxfun(@times, ...
                               exp(1i* (TrajectoryTP.tht + TrajectoryTP.phi)),...
                               geom.L2);

  % initialize the output matrix as NaN's
  distanceTENSOR = nan([nCobra nCobra nSteps]); % MxMxN matrix
  collisionTypeTENSOR = zeros(size(distanceTENSOR));
  
  FIB1 = Trajectory(row,:);
  FIB2 = Trajectory(col,:);
  ELB2 = Elbows(col,:);
  
  [distPt2Line collisionType] = pt2linesegment(FIB1, ELB2, FIB2);
  distanceTENSOR(rcvector) = distPt2Line;
  collisionTypeTENSOR(rcvector) = collisionType;
  
  collisionDetected = distanceTENSOR;
  collisionDetected(distanceTENSOR >  geom.minDist) = false;
  collisionDetected(distanceTENSOR <= geom.minDist) = true;
  collisionDetected(isnan(distanceTENSOR)) = false;
  
  % clear out collision type for cases of no collision.
  collisionTypeTENSOR(~collisionDetected) = 0;
  
  output.detected = collisionDetected;
  output.type     = collisionTypeTENSOR;

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
    
  