function plotsetup(type)
% PLOTSETUP: set up plot for poster or default setup
% 'paper','poster','talk','default'

% taken from peterm@ein:~/matlab/Data/20050228_for_LPSC/makeplot.m
if strcmp(type,'paper')
  set(0,'DefaultAxesFontSize',8) % and legend
  set(0,'DefaultTextFontSize',8)
  set(0,'DefaultAxesPosition',[0.12 0.14 0.85 0.84])
  set(0,'DefaultAxesLineWidth',.5)
  set(0,'DefaultLineLineWidth',.5)
  set(0,'DefaultLineMarkerSize',4)
  clf;
end;

if strcmp(type,'poster')
  set(0,'DefaultAxesFontSize',28)
  set(0,'DefaultTextFontSize',28)
  set(0,'DefaultAxesPosition',[.19 .2 .75 .65])
  set(0,'DefaultAxesLineWidth',1)
  set(0,'DefaultLineLineWidth',1.25)
  set(0,'DefaultLineMarkerSize',10)
  clf;
end;

if strcmp(type,'talk')
  set(0,'DefaultAxesFontSize',18)
  set(0,'DefaultTextFontSize',18)
  set(0,'DefaultAxesPosition',[.15 .15 .8 .7])
  set(0,'DefaultAxesLineWidth',1)
  set(0,'DefaultLineLineWidth',1.25)
  set(0,'DefaultLineMarkerSize',10)
  clf;
end;

if strcmp(type,'default')
  set(0,'DefaultAxesFontSize','factory')
  set(0,'DefaultTextFontSize','factory')
  set(0,'DefaultAxesPosition','factory')
  set(0,'DefaultAxesLineWidth','factory')
  set(0,'DefaultLineLineWidth','factory')
  set(0,'DefaultLineMarkerSize','factory')
  clf;
end

