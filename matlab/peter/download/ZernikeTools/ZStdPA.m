function v = ZStdPA(j, p, A, r0)

% Standard zernike in polar coordinate

% Input:
%   (p, A): polar coordinate of a matrix of points
%   j: the Zernike term
%   r0: radius, default to 1
% Output:
%   v: a matrix of the jth Zernike values at given points (p, A) 

if j>231
    warning('Calculation accuracy may not be guaranteed due to limitations of FACTORIAL!');
end

if(nargin==3)
    r0 = 1;
end

[n m q] = FindNMQ(j); % find sub-indices for 

% angular dependence
if q==1
    Ang = sqrt(2)*cos(m*A);
elseif q==-1
    Ang = sqrt(2)*sin(m*A);
else
    Ang = ones(size(A));
end

% radial dependence
Rad = Radial(n, m, p);
Rad(find(p>r0)) = 0; 

v = sqrt(n+1)*Rad.*Ang;
% figure; imagesc(v); colorbar;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Rad = Radial(n, m, p)
% Radial calculates the radial dependence of the Zernike standard
% polynomials

Rad = zeros(size(p));
for s=0:(n-m)/2
    % calculate
    term = factorial(n-s)/factorial(s)/factorial((n+m)/2-s)/factorial((n-m)/2-s)*p.^(n-2*s);
    if mod(s, 2)==1
        Rad = Rad - term;
    else
        Rad = Rad + term;
    end
end




