function Mout = rotate_matrix(Min, axis, angle, units)
% rotate a 3xn matrix of R3 vectors
% usage: Mout = mrotate(Min, axis, angle)
% axis is the rotation axis
% angle is the RH rotation angle in radians

if ~exist('units','var')
  units = 'radians';
end

if (strcmpi(units,'deg'))
  angle = angle * pi/180;
end

if size(axis,1) < size(axis,2)
  axis = axis';
end
axis = axis/(sqrt(dot(axis,axis))); %normalize the axis
cosphi = cos(angle);
sinphi = sin(angle);
axM = zeros(size(Min));
for jj=1:size(Min,2)
  axM(:,jj) = cross(axis,Min(:,jj));
end

Mout = Min * cosphi + axis * axis' * Min * (1-cosphi) + axM * sinphi;
