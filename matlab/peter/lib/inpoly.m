% inpoly.m
% descended from runs.m by  Z.V. Kovarik
% usage: rr = runs(V,P)
% inputs are vectors of complex numbers
% V is the vertices of the polygon
% P are the points of interest
% output rr is a logical array with 
% 1 for points inside the polygon and 
% 0 for points outside the polygon
%
% If you really need speed, use inpoly by Darren Engwirda on matlab
% central.  It uses a line-crossing algorithm that is 10x faster.
% This one takes a hit when the polygon has many sides.


function rr=inpoly(V,P);

V = V(:); %force into a column vector
P  = transpose(P(:)); % force into a row vector
if (V(1) ~= V(end)) % close path if necessary
  V(end+1) = V(1);
end
zz = bsxfun(@minus, V, P);
rr = round(sum(mod(diff(angle(zz)/pi) + 1, 2) - 1) / 2);
rr(find(isnan(rr))) = 0;
rr = logical(rr);


%from google search: polygon interior algorithm matlab
%
%http://groups.google.com/group/sci.math.num-analysis/browse_thread/thread/10bc010d01dd9760/238a61483a25d7a7%23238a61483a25d7a7?sa=X&oi=groupsr&start=0&num=3
%
%from google search: polygon interior algorithm 
%
%http://astronomy.swin.edu.au/~pbourke/geometry/insidepoly/

