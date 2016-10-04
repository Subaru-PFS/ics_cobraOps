function showMovement(trajectories,geom, coll, pos1, pos2)
% demonstrate x,y to theta, phi transformation for a cobra
% phi is the phi-arm angle wrt the theta arm
% r_err is the radius of the error circle, also in arm-lengths!
% (.0021 == 5 microns)
  psym = 'r';
  thick = 4; % line thickness
  rfib = 1.0; % fiber arm, elbow radius

  % indices in axis 1 of coll.detected corresponding to pos1 and pos2
  rcindx = find((coll.row == pos1 & coll.col == pos2) | ...
                (coll.row == pos2 & coll.col == pos1));
  collIndx = find(sum(coll.detected(rcindx,:)));
  
  clf;
  %% plot collision type vs. step
  subplot(414)
  for jj = 1:length(rcindx)
      plot(coll.type(rcindx(jj),:),'o:'); hold on;
  end
  hold off;
  xlabel('step #');
  ylabel('collision type');
  title('collision of fiber to: (1) elbow (2) arm (3) fiber');
  drawnow;
  
  subplot(4,1,[1 2 3]);
  
  for jj = 1:size(trajectories,2)
    % elbow(1) and fiber(2) xy positions
    xy(1) = trajectories(pos1, jj);
    xy(2) = trajectories(pos2, jj);

    trajectoriesTP = XY2TP(trajectories(:,jj) - geom.center, geom.L1, geom.L2);

    elb1XYloc = geom.center(pos1) + geom.L1(pos1) .* exp(1i*trajectoriesTP.tht(pos1));
    elb2XYloc = geom.center(pos2) + geom.L1(pos2) .* exp(1i*trajectoriesTP.tht(pos2));
    
    % elbow and fiber circles in xy space
    % note: defining circles here is much faster than repeated
    % calls to PLOTCIRCLE.
    fiber1tipXY = xy(1) + rfib * exp(1i*2*pi*(0:.01:1));
    fiber1elbXY = elb1XYloc + rfib * exp(1i*2*pi*(0:.01:1));
    fiber2tipXY = xy(2) + rfib * exp(1i*2*pi*(0:.01:1));
    fiber2elbXY = elb2XYloc + rfib * exp(1i*2*pi*(0:.01:1));
    
    if jj == 1
      cmplx(@plotcircle, geom.center(pos2), geom.L1(pos2) + geom.L2(pos2),'g');
      cmplx(@plotcircle, geom.center(pos1), geom.L1(pos1) + geom.L2(pos1),'b');
      xlabel('X [mm]')
      ylabel('Y [mm]')
      title(sprintf('PFI bench coordinates (X,Y) cobras %d and %d',pos1,pos2))
      hold on;
      axis equal;
    end

    if isempty(find(collIndx == jj))
      plot(fiber1tipXY,'b');
      plot(fiber1elbXY,'b--'); 
      plot(fiber2tipXY,'g');
      plot(fiber2elbXY,'g--'); 
    else
      plot(fiber1tipXY,'r','linewidth',thick);
      plot(fiber2tipXY,'r','linewidth',thick);
      plot(fiber2elbXY,'r:','linewidth',thick); 
    end
    drawnow;
    %disp(jj);
  end
  plot(fiber1tipXY,'k','linewidth' ,thick/2);
  plot(fiber2tipXY,'k','linewidth' ,thick/2);
  plot(fiber1elbXY,'k:','linewidth',thick/2);
  plot(fiber2elbXY,'k:','linewidth',thick/2);
  hold off;
