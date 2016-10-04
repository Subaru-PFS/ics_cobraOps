function [dataout, info] = sigma_reject(datain,Nsigma,max_iterations,verbosity)
% Performs iterative sigma rejection on data
% USAGE: [DataOut, Info] = sigma_reject(DataIn, NSigma, Max_Iterations,verbosity)
% INPUT: DataIn is the data vector
%        NSigma is the rejection criteria
%        Max_Iterations (optional, default=10)
%        verbosity > 1 gives output (not implemented)
% OUTPUT: DataOut = result
%         Info.data = index locations of data
%         Info.index = index locations of data
%         Info.good = logical locations of data
%         Info.iterations = number of iterations taken
%         Info.NSigma = rejection criteria specified
%         Info.orig = DataIn

%DEFAULTS
if ~exist('verbosity')
  verbosity = 0;
end
if ~exist('max_iterations')
  max_iterations = 10;
end

GoodData = logical(ones(size(datain)));
index = cumsum(ones(size(datain)));
last_n_rejected = length(datain(GoodData == 0));

for jj=1:max_iterations
  meanval = mean(datain(GoodData));
  sigma   = std(datain(GoodData));
  AbsDiff = abs(datain - meanval);
  NewGoodData = AbsDiff < Nsigma*sigma;
  new_n_rejected = length(datain(NewGoodData == 0));
  % some debugging crap
  if (verbosity >= 2)
    disp(['j = ',num2str(jj)]);
    disp(['mean = ',num2str(meanval)]);
    disp(['mode = ',num2str(mode(datain(NewGoodData)))]);
    disp(['f_rej = ',num2str(new_n_rejected/length(datain) ) ] );
    disp(sprintf('mean +- sigma = %f +- %f',meanval,sigma));
  end
  % if there are no further rejections, then we are done
  if (new_n_rejected == last_n_rejected)
    jj = jj - 1;
    break
  % it doesn't handle a single value well, so this takes care of
  % that case.  rejecting everything is not an option
  elseif ( new_n_rejected == length(datain) )
    break
  % something was rejected, keep on going...
  else
    last_n_rejected = new_n_rejected;
    GoodData = NewGoodData;
  end
end
switch jj
 case max_iterations
  disp('Warning: maximum number of iterations reached');
 case 0
  if (verbosity >= 1) 
    disp('Note: no data were rejected');
  end
end

dataout = datain(GoodData);
info.data = datain(GoodData);
info.index = index(GoodData);
info.mean = meanval;
info.sigma = sigma;
info.good = GoodData;
info.iterations = jj;
info.NSigma = Nsigma;
info.orig = datain;
