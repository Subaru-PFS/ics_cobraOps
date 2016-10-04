function [n m q] = FindNMQ(j)
% given an index of a Zernike, find its n, m and q. When q is 1, the angle
% dependence is cosine. When it is -1, it's sine.

n = ceil((sqrt(8*j+1)-3)/2);
k = j - n*(n+1)/2;
p = mod(n, 2);
if (p == 0)
    m = 2*floor(k/2);
else
    m = 2*floor((k-1)/2)+1;
end

if(m==0)
    q = 0;
elseif (mod(j, 2))
    q = -1; % sin
else
    q = 1; % cos
end
