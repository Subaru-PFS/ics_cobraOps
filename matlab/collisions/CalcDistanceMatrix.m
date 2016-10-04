function output = CalcDistanceMatrix(Targets, geom)%  Centers,  L1, L2,NNmap)
% given Centers (vector) and a set of (targets|fiberpositions), return the
% pt to line distance for each adjacent pair

  [row col] = find(geom.nnMap);
  rcvector  = find(geom.nnMap);
  
  % compute the elbow position for each target
  targets = XY2TP(Targets - geom.center, geom.L1, geom.L2);
  Elbows = Targets - geom.L2 .* exp(1i* (targets.tht + targets.phi));

  % initialize the output matrix as NaN's
  distance = zeros(size(geom.nnMap));
  
  FIB1 = Targets(row);
  FIB2 = Targets(col);
  ELB2 = Elbows(col);
  
  distPt2Line = pt2linesegment(FIB1, ELB2, FIB2);
  distance(rcvector) = distPt2Line;

% $$$   output.dst=distPt2Line;
% $$$   output.rc = [row col];
  output = sparse(distance);
  
  %% verification plots
  if false
    clf;
    subplot(121)
    cmplx(@plotcircle,Centers,L1 + L2,'r');
    hold on;
    plot(Targets,'o');
    plot([Elbows Targets].');
    hold off;
    
    subplot(122)
    dmap = distance;
    dmap(distance > 2) = 1;
    dmap(distance < 2) = 2;
    dmap(isnan(dmap)) = 0;
    
    imagesc(dmap);
    colorbar;
    caxis([0,2])
    axis equal
    
  end
    