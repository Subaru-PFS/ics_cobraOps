%*****************************************************************************************************
%NAME: errorbars.m
%AUTHOR: Peter Mao, based on Andri M. Gretarsson's work
%DATE: 04/02/2013
%
%SYNTAX: errorbars(@PLOT_FUNCTION,X,DX,Y,DY,'colour')
%
%This function acts just like the built-in plotting functions but plots error-bars. The error-bars
%are plotted in the colour given by 'colour'.  'X','DX','Y', and 'DY' are vectors. 'colour' is a 
%one-letter string which must be one of the letters allowed in the built-in matlab function 'plot', 
%to specify the plot colour.  The function exits with "hold" set to the incoming hold state.
%
%This function does not print points in addition to the error bars.  Where the error bars cross, is
%the coordinate point.  This means that if both error bars are exceedingly small complared to the
%coordinate values, the mark will be correspondingly small.  In such situations, it may be better to
%use 'loglog' directly and specify that the error is smaller than the size of the mark.
%
%EXAMPLE:
%
%X  = [1.0 2.0];
%DX = [0.2 0.2];
%Y  = [1.0 2.0];
%DY = [0.25 0.25];
%errorbars(@plot,X,DX,Y,DY,'g')
%
%plots a green cross of width 0.2 and height 0.25 at coordinate (1.0,1.0), and a cross of width 0.2
%and height 0.25 at coordinate( (2.0,2.0), on a linear scale
%
%LAST MODIFIED:  04/03/2013
%*****************************************************************************************************

function handle=errorbars(PLOTFUNCTION,xx,dx,yy,dy,colourstring)

if (nargin < 5)
  disp('usage: errorbars(@PLOTFUNCTION,X,DX,Y,DY,[''color'']');
  disp('ERRORBARS requires at least 5 inputs!');
  return;
end
if ~exist('colourstring'), colourstring = 'b'; end
if (length(colourstring) > 1) 
  disp('ERRORBARS: color string  argument must be a single color character');
  return;
end
if isempty(regexp(colourstring,'[rgbcmykw]'))
  disp('ERRORBARS: valid color arguments are r,g,b,c,m,y,k,w');
  return;
end
XLOG = false;
YLOG = false;
switch func2str(PLOTFUNCTION)
  case 'loglog'
    XLOG = true;
    YLOG = true;
  case 'semilogx'
    XLOG = true;
  case 'semilogy'
    YLOG = true;
end
init_nohold = ~ishold;

%% make all vectors into row vectors:
xx = transpose(xx(:));
dx = transpose(dx(:));
yy = transpose(yy(:));
dy = transpose(dy(:));

%% for log scaled axes, cut out the negative data
xok = xx < Inf;
yok = yy < Inf;
if XLOG
  xok  = xx > 0;
  if (length(find(xok)) < length(xx)), warning('Negative X data ignored'); end;
end
if YLOG
  yok  = yy > 0;
  if (length(find(yok)) < length(yy)), warning('Negative Y data ignored'); end;
end
ok  = xok & yok;

xx = xx(ok); dx = dx(ok);
yy = yy(ok); dy = dy(ok);

xlo = xx - dx;
xhi = xx + dx;
ylo = yy - dy;
yhi = yy + dy;
if XLOG
  xlo_ok = xlo > 0;
else
  xlo_ok = xlo < Inf;
end
if YLOG
  ylo_ok = ylo > 0;
else
  ylo_ok = ylo < Inf;
end
xlo_bad = find(~xlo_ok); % "bad" ones are where the error range goes negative.
ylo_bad = find(~ylo_ok);

% plot all the well-behaved ones 
h = PLOTFUNCTION([xlo(xlo_ok); xhi(xlo_ok)], [yy(xlo_ok); yy(xlo_ok)], colourstring); hold on;
PLOTFUNCTION([xx(ylo_ok); xx(ylo_ok)], [ylo(ylo_ok); yhi(ylo_ok)], colourstring);

if (XLOG | YLOG)
  % high sides of bad ones
  PLOTFUNCTION([xx(xlo_bad); xhi(xlo_bad)], [yy(xlo_bad); yy(xlo_bad)], colourstring);
  PLOTFUNCTION([xx(ylo_bad); xx(ylo_bad)], [yy(ylo_bad); yhi(ylo_bad)], colourstring);

% $$$ % plot neg ranges in error as positive: this is only useful if the ranges are asymmetric, which is not implemented here.
% $$$ plot_neg_range_as_pos = 1;
% $$$ if (plot_neg_range_as_pos)
% $$$   PLOTFUNCTION(-xlo(xlo_bad), yy(xlo_bad), 'r<');
% $$$   PLOTFUNCTION(xx(ylo_bad), -ylo(ylo_bad), 'rv');
% $$$ end

  xlimits = xlim;
  ylimits = ylim;

  % ...and then the bad ones as dotted lines
  PLOTFUNCTION([xlimits(1)*ones(size(xlo_bad)); xx(xlo_bad)], [yy(xlo_bad); yy(xlo_bad)],...
               colourstring,'linestyle',':');
  PLOTFUNCTION([xx(ylo_bad); xx(ylo_bad)], [ylimits(1)*ones(size(ylo_bad)); yy(ylo_bad)] ,...
               colourstring,'linestyle',':');
end

% return to the original hold state.
if init_nohold
  hold off;
end

handle = h(1);