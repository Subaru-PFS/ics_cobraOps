function output=pt2line3(L1,L2,PT)
% calculate distance from PT to line defined by point L1 and L2
% in three dimensions
% see http://mathworld.wolfram.com/Point-LineDistance3-Dimensional.html

vsize = size(PT);
L1 = L1(:);
L2 = L2(:);
PT = PT(:);

DL = L2 - L1;

alpha = (PT - L1)'*DL/(DL'*DL);

nearest_point  = reshape(L1 + alpha*DL,vsize);
dvec = reshape(PT,vsize) - nearest_point;
dist = norm(dvec);
hdist = norm(dvec(1:2));

output.dist    = dist;
output.xydist  = hdist;
if (length(dvec) == 3) 
  output.zdist   = dvec(3);
end
output.v       = dvec;
output.nearest = nearest_point;
