function extract_stationary_fibers(logfile)

if isunix
  system(sprintf('grep Assign_Centroids %s > centroids.txt'));
  system('perl -pi -e ''s/.*?(\d+).*?([\d\.]+).*?([\d\.]+)\]/$1 $2 $3/'' centroids.txt');
else
  warning(['centroids.txt must be defined. use grep & perl to generate ' ...
           'the file']);
end
  
data = loadxyz('centroids.txt','j','x','y');

xy = hist3(simple(data.c),'edges',{1:2048, 1:2048});

%% threshold in hist for identifying a fiber
thresh = max(reshape(xy,1,[]))/2;

[xcenter ycenter] = find(xy > thresh);
fcenter = xcenter + i*ycenter;


for jj=1:length(xcenter);
  
  fib{jj} = find(abs(data.c - fcenter(jj)) < 3);
  thisdata = data.c(fib{jj});
  threesigma = std(thisdata) / sqrt(2) * 3;
  %  fprintf(1,'%02d  %f\n', jj, threesigma);

  figure(jj);
  plot(data.c(fib{jj}),'.');hold on;
  cmplx(@plotcircle,mean(thisdata),threesigma,'r');  % 3-sigma radius circle
  cmplx(@plotcircle,mean(thisdata),0.025,'m'); % 0.025 pixel radius reference circle
  plot(mean(thisdata),'ro','MarkerSize',10); % this is to help find your fiber when you're zoomed out

  axis equal
  ax1 = axis;
  plot(data.c,'.'); % plot everything
  axis(ax1); % but use the range from the fiber of interest.
  legend('fiber position','3\sigma ref.','r = 0.025 pix','location','NW');
  
  hgsave(jj,sprintf('stationaryFiber%d.fig',jj));
end
  
