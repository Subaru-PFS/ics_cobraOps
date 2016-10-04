function value = Value_ZStdPA(p, A, zn, r0)

% Input:
%   (p, A): polar coordinate of a point
%   zn: Nterm zernike coefficients vector
%   r0: define valid data region p<r0
if nargin==3
    r0 = 1;
end
value = zeros(size(p));

Nterm = length(zn);
for i=1:Nterm
    v = zn(i)*ZStdPA(i, p, A, r0);
    value = value + v;
end
