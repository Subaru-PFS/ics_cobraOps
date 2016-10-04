function yzoom(direction)
%zooms the y axis by a factor of sqrt(10)
%usage: yzoom -- zooms in
%       yzoom(-1) zooms out

if (~exist('direction','var') || (direction >=0))
  direction = 1/sqrt(10);
else
  direction = sqrt(10);
end

ylimits = ylim;
ymean = mean(ylimits);
yextent = abs(diff(ylimits))/2;
yextent = yextent * direction;
ylim([ymean - yextent, ymean + yextent]);
