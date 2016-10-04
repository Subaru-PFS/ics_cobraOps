function result = rebin(data_in,XorY,npts)
% rebins data either by averaging ('x'), summing ('y'),
% or summing, then averaging ('m')
% re-bins a vector by N:1
% USAGE: vector_out = rebin(vector_in,['x'|'y'|'m'],N)
% x vectors are averages of the binned components
% y vectors are the sums of the binned components
% m vector is y/npts -- useful when mean of y is more interesting
% than the sum of y
% example:  x.rebin = rebin(x.raw,'x',3);
% warning: this calculates the values of left-over bins

  if ( ~exist('XorY') || ~ischar('XorY') ...
       || ~regexp(XorY,'[xym]') )
    disp('the second argument should be ''x'' or ''y''');
    return;
  end

  n_in  = length(data_in);
  
%  n_out = floor(n_in / npts);
  n_out = ceil(n_in / npts);
  data_out = zeros(1,n_out);
  
  for (jj = 1:n_out-1)
    range = (jj-1)*npts + 1 : jj*npts;
    if (XorY == 'x')
      data_out(jj) = mean(data_in(range));
    elseif (regexp(XorY,'[ym]'))
      data_out(jj) = sum(data_in(range));
    end
  end
  % handle the last bin carefully, in case it is no a full rebin
  remainder = mod(n_in,npts);
  bin_fraction = (n_in - (n_out - 1) * npts) / npts;
  range = (n_out - 1) * npts + 1 : n_in;
  if (XorY == 'x')
    if (remainder == 0)
      data_out(n_out) = mean(data_in(range));
    elseif (remainder > 1)
      % (inter/extra)polates the last bin location from the existing values
      lo = min(data_in(range));
      hi = max(data_in(range));
      step = (hi - lo) / (max(range) - min(range));
      data_out(n_out) = (remainder * step)/bin_fraction*0.5 + lo - step/2;
    elseif (remainder == 1)
      % extrapolates from previous bins
      if (n_out > 2)
	data_out(n_out) = 2 * data_out(n_out-1) - data_out(n_out-2);
      else
	data_out(n_out) = NaN;
      end
    end
  elseif (regexp(XorY,'[ym]'))
    data_out(n_out) = sum(data_in(range)) / bin_fraction;
  end
  if (XorY == 'm')
    data_out = data_out/npts;
  end
    
  result = data_out;
  
  return