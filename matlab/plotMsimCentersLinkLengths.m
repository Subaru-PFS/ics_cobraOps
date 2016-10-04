function plotMsimCentersLinkLengths(data, pidORpnum)
% plot function to go with processMsimCentersLinkLengths_no_plots
% data is the output of processMsimCentersLinkLengths_no_plots
% when 2nd argument is negative, it is an array position (pnum)
% when positive, it is a pid.

%% get array number from pid

  pids = [data.tht.pid];
  
  if pidORpnum < 0
    pnum = -pidORpnum;
  else
    pnum = find(pids == pidORpnum);
  end
  %% THETA

    ZZ = data.tht(pnum).centroids - data.tht(pnum).center;
    angles = angle(ZZ);
    radii  = abs(ZZ);
    ok = data.tht(pnum).ok;

    figure(10)
    axis([0 2000 0 2000]);
    daspect([1,1,1]);
    hold on;
    plot([data.tht.center],'rx')
    plot([data.tht.centroids], 'bo');
    plot(data.tht(pnum).centroids,'gx');
    plot([data.phi.center],'rx')
    plot([data.phi.centroids], 'bo');
    plot(data.phi(pnum).centroids,'gx');
    hold off;
    xlabel('X [pix]');
    ylabel('Y [pix]');
    title(sprintf('\\Theta and \\Phi(pId %d)',data.tht(pnum).pid));

    figure(11) 
    plot(angles/(2*pi),radii, 'bx');
    hold on;
    plot(angles(ok)/(2*pi), radii(ok),'go');%,'MarkerFace','g');
    hold off;
    xlabel('angle/\tau');
    ylabel('radius [pix]');
    title(sprintf('\\Theta PID %d',data.tht(pnum).pid));

    figure(12)
    hist(radii(ok),50);
    xlabel('radius [pix]');
    title(sprintf('\\Theta PID %d',data.tht(pnum).pid));
    
    %% PHI    
    ZZ = data.phi(pnum).centroids - data.phi(pnum).center;
    angles = angle(ZZ);
    radii  = abs(ZZ);
    ok = data.phi(pnum).ok;

% $$$     figure(20)
% $$$     axis([0 2000 0 2000]);
% $$$     daspect([1,1,1]);
% $$$     hold on;
% $$$     plot([data.phi.center],'rx')
% $$$     plot([data.phi.centroids], 'bo');
% $$$     plot(data.phi(pnum).centroids,'gx');
% $$$     hold off;
% $$$     xlabel('X [pix]');
% $$$     ylabel('Y [pix]');
% $$$     title(sprintf('\\Phi (pId %d)',data.phi(pnum).pid));
    
    figure(21) 
    plot(angles/(2*pi),radii, 'bx');
    hold on;
    plot(angles(ok)/(2*pi), radii(ok),'go');%,'MarkerFace','g');
    hold off;
    xlabel('angle/\tau');
    ylabel('radius [pix]');
    title(sprintf('\\Phi PID %d',data.phi(pnum).pid));

    figure(22)
    hist(radii(ok),50);
    xlabel('radius [pix]');
    title(sprintf('\\Phi Pid %d',data.tht(pnum).pid));
    return;
    
