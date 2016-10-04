function [value]=ZStdXY(nZern,x,y)

p = sqrt(x.*x + y.*y);
A = atan2(y, x);

value = ZStdPA(nZern, p, A);