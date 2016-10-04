function PlotHandle = refline(x0,y0,m,linetype)
% REFLINE.M plots a reference line on the current figure
%
% USAGE: [plot_handle =] refline(x,y,slope,linetype)
%
% example: refline(1,1,.5,'g--')
% use slope = (NaN|Inf) for a vertical line

  if isnan(m) || isinf(m)
    xxplot = [x0 x0];
    yyplot = ylim;
  elseif (m == 0)
    xxplot = xlim;
    yyplot = [y0 y0];
  else
    xx = xlim;
    yy = ylim;
    if (m < 0) % must keep y sorted properly
      yy = fliplr(yy);
    end
    
    dx = xx - x0;
    dy = yy - y0;
    
    yycalc = dx * m + y0;
    xxcalc = dy / m + x0;
    
    %  valid points are: 
    %  xx(1) yycalc(1) ||  xxcalc(1) yy(1)
    %  xx(2) yycalc(2) ||  xxcalc(2) yy(2)
    
    if (xxcalc(1) < xx(1))
      xxplot(1) = xx(1);
      yyplot(1) = yycalc(1);
    else
      xxplot(1) = xxcalc(1);
      yyplot(1) = yy(1);
    end
    if (xxcalc(2) > xx(2))
      xxplot(2) = xx(2);
      yyplot(2) = yycalc(2);
    else
      xxplot(2) = xxcalc(2);
      yyplot(2) = yy(2);
    end
  end
  do_not_hold = ~ishold;
  hold on;
  PlotHandle = plot(xxplot,yyplot,linetype);
  if (do_not_hold)
    hold off;
  end
  