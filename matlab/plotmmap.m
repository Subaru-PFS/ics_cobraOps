function plotmmap(plotName, motdir, newMap, oldMap, noNegMoves, subfig)
plothandle = subplot(2,2,subfig);
title(plotName);
hold on;
startAngleH = motdir.startAngle;
finishAngleH = motdir.finishAngle;
mmapH = motdir.mmap;

secondHandle = plot(startAngleH,mmapH,'*r'); % Move Start Angles
plot(finishAngleH,mmapH,'*r'); % Move Finsh Angles
labels = cellstr([num2str(motdir.targetNo-1, '%d_'), num2str(motdir.iteration)]);
%keyboard
if(~noNegMoves)
text(finishAngleH,mmapH, labels, 'VerticalAlignment','bottom', ...
                               'HorizontalAlignment','right');
end
thirdHandle = plot([startAngleH, finishAngleH]',[mmapH,mmapH]','r'); % Just the line

h2 = plot(oldMap(1,:)*180/pi, oldMap(2,:), 'go-', 'linewidth', 3);  % old map
h1 = plot(newMap(1,:), newMap(2,:),'bo-', 'linewidth',3); % new map
h3 = plot(newMap(1,:), newMap(2,:) + newMap(3,:), 'bo-');
plot(newMap(1,:), newMap(2,:) - newMap(3,:), 'bo-');

leg{1} = 'new map';
leg{2} = 'old map';
legend([h1,h2],leg{:});
ylabel('move size in [deg/step]')
if(subfig <3)
    if(noNegMoves)
    axis([0 360 0.0 0.3]);
    else
    axis([0 360 -1.0 0.3]);
    end
    xlabel('theta in [deg]');
else
    if(noNegMoves)
        axis([0 200 0.0 0.3]);
    else
        axis([0 200 -1.0 0.3]);
    end
    xlabel('phi in [deg]');
end

set(findall(plothandle,'type','text'),'FontSize',14,'fontWeight','bold');
set(findall(plothandle,'type','axes'),'FontSize',14,'fontWeight','bold');

end

% dataX = rand(1000,1);
% dataY = rand(1000,1);
% cmap = jet(length(1000));
% figure(12)
% scatter(dataX, dataY, 10, cmap, 'filled')
% 
% 
% x = rand(1000,1);
% y = rand(1000,1);
% distances = ((x-.5).^2 + (y-0.5).^2).^0.5;
% % Normalize - divide my sqrt(maxX^2 + maxY^2)
% distances = distances / sqrt(.5^2 + .5^2);
% [sortedDistances sortIndexes] = sort(distances);
% % Arrange the data so that points close to the center
% % use the blue end of the colormap, and points 
% % close to the edge use the red end of the colormap.
% xs = x(sortIndexes);
% ys = y(sortIndexes);
% cmap = jet(length(x)); % Make 1000 colors.
% scatter(xs, ys, 10, cmap, 'filled')
% grid on;
% title('Points where color depends on distance from center', ...
% 	'FontSize', 30);
% % Enlarge figure to full screen.
% set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
% % Give a name to the title bar.
% set(gcf,'name','Demo by ImageAnalyst','numbertitle','off') 