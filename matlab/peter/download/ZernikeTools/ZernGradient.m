function [gradx grady] = ZernGradient(zerncoeff, x, y)
% Input:
%   zerncoeff: a COLUMN vector containing coefficients for Zernike standard
%   polynomials
%   (x, y): matrices containing normalized (x, y) coordinate of points
%   where gradients are to be calculated


numz = length(zerncoeff);

[CMatX CMatY] = GenerateCoeffMatrix(numz, numz);
Zrep = repmat(zerncoeff, 1, numz);
C_dPhi_dx = sum(Zrep.*CMatX); % row vector of zernike coefficients for x slope
C_dPhi_dy = sum(Zrep.*CMatY); % row vector of zernike coefficients for y slope

p = sqrt(x.*x + y.*y);
A = atan2(y, x);
gradx = Value_ZStdPA(p, A, C_dPhi_dx);
grady = Value_ZStdPA(p, A, C_dPhi_dy);
figure; imagesc(gradx); title('X component of the gradient'); axis equal; colorbar;
figure; imagesc(grady); title('Y component of the gradient'); axis equal; colorbar;


