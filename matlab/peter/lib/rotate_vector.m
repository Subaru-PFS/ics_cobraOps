function vout = rotate_vector(vin, axis, angle, units)
% rotate vector
% vin is the operand
% axis is the direction of the rotation axis
% angle is the right-handed rotation in radians

if ~exist('units','var')
  units = 'radians';
end

if (strcmpi(units,'deg'))
  angle = angle * pi/180;
end

axis = axis/(sqrt(dot(axis,axis))); %normalize the axis
cosphi = cos(angle);
sinphi = sin(angle);
vout = vin * cosphi + axis * dot(axis,vin) * (1-cosphi) + ...
    cross(axis,vin) * sinphi;
