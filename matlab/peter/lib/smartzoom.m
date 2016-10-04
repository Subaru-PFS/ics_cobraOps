function smartzoom(percentage, pos_percentage)
% SMARTZOOM: Zooms away the peaks in the data.  Sort of like "axis
% auto".  Good for looking at noise.
%
% usage: >> smartzoom     % run with default settings
%        >> smartzoom(1)  % zoom view to reject 1% of the data, .5% lost
%        on either end
%        >> smartzoom(1,2) % zoom view to reject 1% low, 2% high.
%
% default: -Zooms to lose 5% symmetrically (2.5% off top and bottom)
%          -stretches the resulting range by 20%, so interesting
%          data is not sitting on the axes.
% 
% smartzoom(0) is tighter than "ylim auto"

if ~exist('percentage','var')
  percentage = 5;
end
if ~exist('pos_percentage','var')
  highcut = 100-percentage/2;
  lowcut  = percentage/2;
else
  highcut = 100-pos_percentage;
  lowcut  = percentage;
end

pstruct = get(gca);

if (pstruct.Children == 1) % single plot is easy
  childOBJ = get(pstruct.Children);
  ydata = childOBJ.YData;
else % multiple data sets require some special handling
  ydata = []; 
  for jj=1:length(pstruct.Children)
    childOBJ = get(pstruct.Children(jj));
    if isfield(childOBJ,'YData')
      ydata = [ydata childOBJ.YData];
    end
  end
end

[NN, XX] = hist(ydata,1000);
DX = XX(2)-XX(1); % the histogram binwidth
dist = cumsum(NN)/sum(NN)*100; % cumulative distribution

xmax = min(find(dist>highcut)); %first point that crosses the high threshold
xmin = max(find(dist<lowcut));  %last point below the low threshold
if isempty(xmax)
  xmax = length(NN);
end
if isempty(xmin)
  xmin = 1;
end
ymin = XX(xmin)-DX/2;  % get ymin/ymax, shift to the bin edges
ymax = XX(xmax)+DX/2;
ymid = (ymin + ymax)/2;
halfrange = (ymax - ymin)/2 * 1.2; %stretch the range by 20%

ylim([ymid - halfrange,  ymid + halfrange]);
