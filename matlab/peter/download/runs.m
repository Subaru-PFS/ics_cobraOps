% runs.m     calculates the index of a point P with respect
%            to a closed polygonal path given by n points.
% ******>    Input:  XY , an n-by-2 matrix XY (the rows
%                         defining the points),
%                    P,   a 2-component vector.
% ******>    Output: r,   the net number of runs of XY around P
%            It is assumed that the last point is connected
%            to the first point.
% ******>    Call: r=runs(XY,P)
%            By Z.V. Kovarik, May 3, 1996
%
function rr=runs(XY,P);
% $$$ i=sqrt(-1);
% $$$ P=P(1)+i*P(2);
% $$$ Z=XY(:,1)+i*XY(:,2);
% $$$ n=length(Z);
% $$$ Z=[Z;Z(1)]-P;             % closing up the path
% $$$ Z=Z(1+(1:n))./Z(1:n);     % complex rotations - unscaled
% $$$ r=round(sum(atan(imag(Z)./(abs(Z)+real(Z))))/pi);

XY = XY(:); %force into a column vector
P  = transpose(P(:)); % force into a row vector
if (XY(1) ~= XY(end)) % close path if necessary
  XY(end+1) = XY(1);
end
zz = bsxfun(@minus, XY, P);
rr = round(sum(mod(diff(angle(zz)/pi) + 1, 2) - 1) / 2);

% end of the function file runs.m

%from google search: polygon interior algorithm matlab
%
%http://groups.google.com/group/sci.math.num-analysis/browse_thread/thread/10bc010d01dd9760/238a61483a25d7a7%23238a61483a25d7a7?sa=X&oi=groupsr&start=0&num=3
%
%from google search: polygon interior algorithm 
%
%http://astronomy.swin.edu.au/~pbourke/geometry/insidepoly/

