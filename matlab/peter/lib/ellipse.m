function output=ellipse(a,b,angle,center,npts)
% generate an ellipse

if ~exist('center','var'), center = 0; end;
if ~exist('npts','var'), npts = 100; end;

tht = 2*pi*(0:npts-1)/npts;

RR = a*b./sqrt((b*cos(tht)).^2 + (a*sin(tht)).^2);

XY = RR .* exp(i*(tht + angle)) + center;

output = XY;