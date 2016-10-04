function [CMatX CMatY] = GenerateCoeffMatrix(numZ, numDZ)

% generate coefficient matrix for Zernike gradient
% this is for calculating slope of a phase function defined by Zernike
% standard phase

% Input:
%       numZ: number of Zernike terms used in expression of Zernike
%       gradient
%       numDz: number of Zernike gradient terms

% xslopeZ = CMatX * Z
% yslopeZ = CMatY * Z

if(nargin == 0)
    numZ = 50;
    numDZ = 50;
end

CMatX = [];
CMatY = [];

for i=1:numDZ
    coeffX = DerivX_Zern(i, numZ);
    CMatX = [CMatX; coeffX'];
    coeffY = DerivY_Zern(i, numZ);
    CMatY = [CMatY; coeffY'];
end
    
    
