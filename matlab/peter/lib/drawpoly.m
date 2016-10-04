function [x,y] = drawpoly()
% DRAWPOLY: draw a polygon on a plot graph
% 
% USAGE: [x,y] = drawpoly;
% 
% use left button to choose points, right button for last point
% function automatically closes the polygon
  
  not_held = ~ishold;
  hold on;
  button = 1;
  n = 0;
  while button == 1
    [xi, yi, button] = ginput(1);
    plot(xi,yi,'ro')
    n = n+1;
    x(n) = xi;
    y(n) = yi;
  end
  x(n+1) = x(1);
  y(n+1) = y(1);
  plot(x,y,'r');
  if not_held
    hold off;
  end