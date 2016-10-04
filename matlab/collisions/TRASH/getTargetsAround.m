function target = getTargetsAround(center)
 
  numtrg = 1000; % Specify number of targets to generate
 
  L1 = 2.375;
  L2 = 2.375;
  
  % set range of phi arm to be ~5deg inside of 0-180.
  KeepOutAngle = 0.1;

  phiOut = KeepOutAngle;
  phiIn  = pi - KeepOutAngle;
  
  % Inner keepout radius
  Rmin = abs(L1 + L2 * exp(i*phiIn));

  % Outer reach
  Rmax = abs(L1 + L2 * exp(i*phiOut));
  
  % angular coordinate of target
  THT = rand(1,numtrg)*2*pi;

  % radial coordinate of target
  RDS = -ones(1,numtrg);
  
  while(true)
    % "indices of radii to assign"
    jj_R_assign = find(RDS < Rmin);
    n_assign = length(jj_R_assign);
    % asssign radii, choosing uniformly over a disk of radius Rmax
    % or end if there's nothing to do.
    if n_assign > 0
      RDS(jj_R_assign) = Rmax * sqrt(rand(1,n_assign));
    else
      break;
    end
  end
  
  target = center + RDS .* exp(1i*THT); 
end