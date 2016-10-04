function page_setup(fig_handle, width, height, orientation, center)
% PAGE_SETUP.M: resizes a plot on A4 paper (20.99 x 29.70 cm)
%
% USAGE: page_setup(gcf, width, height, orientation, center)
% width and height are that of a figure
% orientation is for paper and {'portrait'} or 'landscape'
% units are in centimeters
% set center to 1 to center plot, otherwise it will show up in
% the upper left corner
%
% 8"x6" is the matlab default print size
%
% example: page_setup
% example: page_setup(gcf)
% example: page_setup(gcf, 8, 10, 'portrait', 1)
% example: page_setup(gcf, 8, 10, 'portrait', 1)
%
% 03/04/2006: plotsize.m was by Peter Mao for US letter
% 12/20/2007: modified by Tak Kunihiro for A4
% 02/18/2008: modified by Tak Kunihiro for Letter


%--------------------------------------------------
% handle arguments
%--------------------------------------------------
%height_my_paper = 29.70;
%width_my_paper  = 20.99;
height_my_paper = 11;
width_my_paper  = 8.5;

if ~exist('fig_handle','var')
  fig_handle = gcf;
end
if ~exist('width','var')
%  width = 8;
  width = width_my_paper;
end
if ~exist('height','var')
%  height = 6;
  height = width / width_my_paper * height_my_paper;
end
if ~exist('orientation','var')
  orientation = 'portrait';
end
if ~exist('center','var')
  center = 0;
end
%--------------------------------------------------



%--------------------------------------------------
% window setup
%     maximize window size with aspect ratio of the plot
%--------------------------------------------------
rect_screen = get(0, 'ScreenSize');
height_win = rect_screen(4) * 0.80;
width_win  = height_win / height * width; % proportional to the plot
left_win   = rect_screen(2)+5;
bottom_win = height_win * 0.05;
set(fig_handle, 'Position', [ left_win bottom_win width_win height_win ]);
%--------------------------------------------------



%--------------------------------------------------
% paper setup
%--------------------------------------------------
if strcmp(orientation,'portrait')
  paperwidth = width_my_paper;
  paperlength= height_my_paper;
elseif strcmp(orientation,'landscape')
  paperwidth = height_my_paper;
  paperlength= width_my_paper;
end

if (center)
  xmargin = (paperwidth  - width )/2;
  ymargin = (paperlength - height)/2; % no /2, want to push it to the
                                  % top of the page
else
  xmargin = 0;
  ymargin = paperlength - height;
end

orient(orientation);
set(fig_handle, 'PaperPositionMode','manual');
%set(gcf, 'PaperUnits', 'inches');
set(fig_handle, 'PaperUnits', 'inches');
set(fig_handle, 'PaperType', 'usletter');
%set(fig_handle, 'PaperUnits', 'centimeters');
%set(fig_handle, 'PaperType', 'A4');
set(fig_handle, 'PaperPosition',[ xmargin ymargin width height ]);
%--------------------------------------------------




return
