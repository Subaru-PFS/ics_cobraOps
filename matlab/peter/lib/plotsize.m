function plotsize(fighandle, width, height, orientation, center)
% PLOTSIZE.M: resizes a plot on US letter paper
%
% USAGE: plotsize(gcf, width, height, orientation, center)
% units are in inches
% orientation is {'portrait'}, 'landscape', or 'electronic'
% "electronic" orientation is papersize = figuresize
% set center to 1 to center plot, otherwise it will show up in
% the lower left corner
%
% 8"x6" is the matlab default print size
%
% example: plotsize(gcf, 8, 10, 'portrait', 1)
%
% plotsize defaults to plotsize(gcf,6,4.5,'electronic',0);

if ~exist('fighandle','var')
  fighandle = gcf;
end
if ~exist('width','var')
  width = 6;%8;
end
if ~exist('height','var')
  height = 4.5;%6;
end
if ~exist('orientation','var')
  orientation = 'electronic';%'portrait';
end
if ~exist('center','var')
  center = 0;
end

if strcmp(orientation,'portrait')
  paperwidth = 8.5;
  paperheight= 11.0;
elseif strcmp(orientation,'landscape')
  paperwidth = 11.0;
  paperheight= 8.5;
elseif strcmp(orientation,'electronic')
  paperwidth = width;
  paperheight = height;
end

if (center)
  xmargin = (paperwidth - width)/2;
  ymargin = (paperheight - height)/2; % no /2, want to push it to the
                                  % top of the page
else
  xmargin = 0;
  ymargin = paperheight - height;
end

set(fighandle, 'PaperPositionMode','manual');
set(gcf, 'PaperUnits', 'inches','PaperSize',[paperwidth paperheight]);
set(gcf,'PaperPosition',[ xmargin ymargin width height ]);
return
