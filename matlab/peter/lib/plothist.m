function plothist(filename,varname)
% plots data into a histogram....
% taken from 20060120_delta18_em/plothist.m

%q=loadlogger(filename);
q = readcdf(filename);
%rawdata=q.time; % have to declare this or eval won't work
%testcmd  = ['rawdata = q.',varname,';']
%eval(testcmd);
rawdata = q.(varname).data;

n_samples = length(rawdata);
nbins = round(10*log10(n_samples));
[qh,xh] = hist(rawdata,nbins);
binwidth = xh(2)-xh(1);

qh_mean = mean(rawdata);
qh_std = std(rawdata);
params(2) = qh_mean;
params(3) = qh_std;

params(1) = binwidth*n_samples;
disp(['using params ', num2str(params)])
qhfit = gaussian(xh,params);

figure(101)
subplot(211)
plot(rawdata)
title(strrep(filename,'_','\_'))

subplot(212)
bar(xh,qh)
hold on
plot(xh,qhfit,'r')
hold off
delta = sum(chisq(qh,qhfit,qh))/(length(xh)-1);
legend('distribution',['\chi^2_\nu = ',num2str(delta)]);
disp(['reduced chi^2 = ',num2str(delta)]);
title(strrep(filename,'_','\_'))
return

function y = gaussian(x,a)
  % a = [amplitude, mean, sigma]
  amp = a(1);
  mean = a(2);
  sigma= a(3);
  
  y = amp/(sigma*sqrt(2*pi)) * exp(-1*(x-mean).^2/(2*sigma^2));
  return

function answer = chisq(y,yfit,variance)

  n_elements = length(y);

  if nargin() < 3
    variance = ones(size(y));
  end
  sn2 = (y-yfit).^2/n_elements;
  answer = n_elements*sn2./max(variance,1);
  return

