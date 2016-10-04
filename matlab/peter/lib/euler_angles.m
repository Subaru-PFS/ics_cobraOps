function Mout = euler_angles(Min, theta, phi, psi, units);
% rotate matrix given body-fixed Euler angles
% Mout = euler_angles(Min, theta, phi, psi, units);

if ~exist('units','var')
  units = 'radians';
end

if (strcmpi(units,'deg'))
  theta = theta * pi/180;
  phi = phi * pi/180;
  psi = psi * pi/180;
end

M1 = rotate_matrix(Min, [0 1 0], theta);
M2 = rotate_matrix(M1, M1(:,1), phi);
Mout = rotate_matrix(M2,M2(:,3),psi);
