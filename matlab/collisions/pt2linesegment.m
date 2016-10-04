function [dist solntype] = pt2linesegment( PT, L1, L2 )
% calculate distance from PT to a line segment defined by point L1 and L2
% in two dimensions
% all coordinates are complex numbers.
  
%% translate line segment and point so that L1 is on the origin
  line1 = L2 - L1;
  point = PT - L1;

  %% rotate line & points so that the line lies on the +X axis
  rotated_point = point .* exp(-i*angle(line1));
  
  %% define 3 regions: left of origin, inside, right of L2 end
  %% point.
  XX = real(rotated_point);
  X2 = abs(line1); % X coordinate of L2 in T+R coordinates
  REG1 = find(XX <= 0);
  REG2 = find(XX > 0 & XX < X2);
  REG3 = find(XX >= X2);
  NANS = find(isnan(XX) | isnan(X2));
  
  dist(REG1) = abs(L1(REG1) - PT(REG1));
  dist(REG2) = abs(imag(rotated_point(REG2)));
  dist(REG3) = abs(L2(REG3) - PT(REG3));
  dist(NANS) = nan(size(NANS));
  dist = reshape(dist, size(PT));
  
  solntype(REG1) = 1;
  solntype(REG2) = 2;
  solntype(REG3) = 3;
  solntype(NANS) = 0;
  solntype = reshape(solntype, size(PT));
  
% $$$   dist(1:numel(rotated_point)) = abs(imag(rotated_point)); 
% $$$   dist(real(rotated_point) < 0) = distabsBefore(real(rotated_point) < 0);
% $$$   dist(real(rotated_point) > abs(line1)) = distabsBehind(real(rotated_point) > abs(line1));
 
end