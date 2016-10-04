function showMovementNN(trajectories,geom, coll, pos1, targets)
% show trajectories for Pos and all neighbors

  cobraPos  = [pos1 find(geom.nnMap(pos1,:))];
  localIndx = find(ismember(coll.row, cobraPos) & ismember(coll.col, cobraPos));
  localColl = coll.detected(localIndx, :);
  collInTrajStep = sum(localColl);

  isMoving = find([true diff(abs(sum(trajectories(cobraPos,:)))) ~= 0]);

  % loop over all trajectory steps where something moves
  for kk = 1:length(isMoving)
      jj = isMoving(kk);
      xy = trajectories(cobraPos,jj);
      trajectories_cobraPos = XY2TP(trajectories(cobraPos,jj) - geom.center(cobraPos),...
                             geom.L1(cobraPos), geom.L2(cobraPos));
      elbows = geom.center(cobraPos) + geom.L1(cobraPos) .* exp(1i*trajectories_cobraPos.tht);

      % phi arm
      set(gca,'ColorOrderIndex',1)
      plot([elbows xy].','o-','linewidth',1);
      if jj == 1
          hold on;
          cmplx(@plotcircle, geom.center(cobraPos), geom.rMax(cobraPos),'k:');
          plot(targets(cobraPos),'mo','MarkerFace','m');
          title(sprintf('Trajectories around positioner %d', pos1));
          xlabel('X [mm]');
          ylabel('Y [mm]');
      end
      % theta arm
      set(gca,'ColorOrderIndex',1)
      plot([geom.center(cobraPos) elbows].','-');
      % patrol regions
      if collInTrajStep(jj) % a collision is detected on this step
          % get the localColl row and col numbers in xy
          kk = find(localColl(:,jj)); 
          rr = find(ismember(cobraPos, coll.row(localIndx(kk))));
          cc = find(ismember(cobraPos, coll.col(localIndx(kk))));
          plot(xy(rr),'ro','MarkerFace','r');
          plot([elbows(cc) xy(cc)].','ro-','MarkerFace','r','linewidth',3);
          cmplx(@plotcircle, xy(rr), geom.rf, 'r');
          cmplx(@plotcircle, xy(cc), geom.rf, 'r');
          cmplx(@plotcircle, elbows(cc), geom.rf, 'r');
      end
      drawnow;
  end
  cmplx(@text, geom.center(cobraPos), cellstr(num2str(cobraPos')));
  hold off;
