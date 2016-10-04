function [value_X, value_Y] = S_xy(index, x, y)
% calculate value of the vector polynomial S with index "index" at point
% (x, y)
% output:
%   value_X: x component of the value
%   value_Y: y component of the value

if index>80 || index<=0
    error('Out of range!');
end

value_X = zeros(size(x));
value_Y = zeros(size(y));

[Z_Sjx C_Sjx Z_Sjy C_Sjy] = FindSj(index);
p = sqrt(x.*x + y.*y);
A = atan2(y, x);

nx = length(Z_Sjx); % number of Zernikes in x-component
for i = 1:nx
    value_X = value_X + C_Sjx(i)*ZStdPA(Z_Sjx(i), p, A);
end

ny = length(Z_Sjy); % number of Zernikes in y-component
for i = 1:ny
    value_Y = value_Y + C_Sjy(i)*ZStdPA(Z_Sjy(i), p, A);
end
