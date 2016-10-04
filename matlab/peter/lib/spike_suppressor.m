function [dataout, info] = spike_suppressor(datain,window,Nsigma,reject_crit,verbosity)
% remove transient spikes from data.  sigma rejects in a window.
% USAGE: [DataOut, Info] = spike_suppressor(DataIn, window, [NSigma,reject_crit,verbosity])
% REQUIRES: SIGMA_REJECT, ALLSTATS
% INPUT: DataIn = data vector
%        window = width of comparison window
%        NSigma = sigma rejection criteria.  Default is (0.5*sqrt(window))
%        reject_crit = 0 rejects only the worst data, 1 rejects
%        everything. default is 0.85.
%        verbosity > 1 for debugging
% OUTPUT: DataOut = result
%         Info.data = result
%         Info.index = index locations of data
%         Info.good = logical goodness (1=keep)
%         info.goodness = real number goodness
%         Info.NSigma = sigma rejection criteria specified
%         Info.orig = original data
%
% SPIKE_SUPPRESSOR runs SIGMA_REJECT on a sliding region (of width
% "window") across the data set, looking for spikes.  A line fit is
% subracted from the region before the rejection algorithm
% operates, so slowly varying data will not be rejected.  Spikes
% should appear in the majority of the comparisons as the analysis
% region moves through the spike.  Each time a data point is
% rejected by SIGMA_REJECT, it's "goodness" (starting at 1) is
% reduced.  If it is always rejected, then its "goodness" will go
% to zero.  "reject_crit" determines the level of "goodness" below
% which we actually reject the data.

%DEFAULTS
if ~exist('verbosity','var')
  verbosity = 0;
end
if ~exist('max_iterations','var')
  max_iterations = window;
end
if ~exist('reject_crit','var')
  reject_crit = .85;
end
if ~exist('Nsigma','var')
  Nsigma = 0.5*sqrt(window);
end

%SANITY CHECKS
quit_early = 0;
if window > length(datain)
  warning('window width > # data points');
  quit_early = 1;
end
if window < 3
  warning('window too narrow');
  quit_early = 1;
end
if Nsigma^2 > window
  warning('Nsigma too large for window width');
  quit_early = 1;
end
if quit_early
  dataout = datain;
  info.data = datain;
  info.index = 1:length(datain);
  info.orig = datain;
  return
end

%SET UP ARRAYS AND VARIABLES
GoodData = ones(size(datain));
index = 1:length(datain);
window1 = window - 1;
iterations = length(datain) - window1;
current_fig = gcf;
% kk is the window location index
for kk=1:iterations
  roi = kk:kk + window1;
  temp = datain(roi);
  % FIT DATA SEGMENT TO A LINE
  P = polyfit(roi,temp,1);
  temp = temp - polyval(P,roi);
  [temp_rej temp_info] = sigma_reject(temp,Nsigma,max_iterations);
  % if a point rejects, then proportionally reduce the goodness of
  % the point
  GoodData(kk:kk+window1) = GoodData(kk:kk+window1) ...
      - (~temp_info.good/min([window iterations]) );  
  % DEBUGGING PLOTS
  if (verbosity > 0)
    stats = allstats(temp);
    figure(99)
    subplot(221)
    plot(roi,temp,'.'); hold on;
    refline(mean(roi),stats.mean+stats.std*Nsigma,0,'.--');
    refline(mean(roi),stats.mean-stats.std*Nsigma,0,'.--');
    hold off;
    ylabel('ROI data - line fit');
    subplot(223)
    plot(roi,temp_info.good);
    ylabel('ROI goodness');
    subplot(224)
    plot(GoodData);
    refline(length(datain)/2, reject_crit, 0, '--');
    ylabel('data quality')
    figure(current_fig)
  end
end

Good = GoodData > reject_crit;

dataout = datain(Good);
info.data = datain(Good);
info.index = index(Good);
info.good = Good;
info.goodness = GoodData;
info.reject_crit = reject_crit;
info.NSigma = Nsigma;
info.orig = datain;
  
return
