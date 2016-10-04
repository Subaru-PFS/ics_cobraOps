function [value_X, value_Y] = T_xy(index, x, y)
% calculate value of the vector polynomial T with index "index" at point
% (x, y)
% output:
%   value_X: x component of the value
%   value_Y: y component of the value

% calculate S
[value_X, value_Y] = S_xy(index, x, y);

% get T
temp = value_X;
value_X = value_Y;
value_Y = -temp;
