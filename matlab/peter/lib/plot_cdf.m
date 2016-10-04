function output=plot_cdf(list_of_errors,psym, gausscmp, gausscmp_psym)
% plot the cumulative distribution of a list of errors
% usage PLOTHANDLE = PLOT_CDF(list_of_errors, [PSYM1, SHOW_GAUSSIAN_logical, PSYM2])

  if ~exist('psym','var'), psym = 'b'; end;

  ISCOMPLEX = false;
  if any(~isreal(list_of_errors))
    ISCOMPLEX = true;
  end
  
  sorted_list = sort(abs(list_of_errors));
  
  nn = length(list_of_errors);
  
  yy = 100 * (1:nn)/nn;
  output = semilogx(sorted_list, yy, psym);
  xlabel('[mm]');
  ylabel('culumative %');
  
  if exist('gausscmp','var')
    if ~exist('gausscmp_psym','var'), gausscmp_psym = 'r--'; end
    hold on;
    if ISCOMPLEX
      sigma = std(list_of_errors) * sqrt(0.5);
      gg = sigma * (randn(nn,1) + i*randn(nn,1));
    else
      sigma = std(list_of_errors);
      gg = sigma * randn(nn,1);
    end
    gs = sort(abs(gg));
    plot(gs,yy,gausscmp_psym);
    hold off;
  end
end