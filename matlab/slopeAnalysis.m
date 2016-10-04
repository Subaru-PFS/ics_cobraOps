clear all
%close all
clc

slopeFiles = dir2cell('slopes_*');
numFiles = length(slopeFiles);

fh1 = figure('Name','Slope Comparison');hold on
fh2 = figure('Name','Normalized Slope Comparison');hold on
fh3 = figure('Name','Slope Deviation from Initial Conditions');hold on

%% Initialize slopeStruct with "slopes*" output file from adjustMaps.m
for ii=1:numFiles
    tempData = dlmread(slopeFiles{ii},',',1,0);
    slopeStruct(ii).file = slopeFiles{ii};
    pId = tempData(:,1);
    slopeStruct(ii).pId = pid;
    for jj=1:length(pId)
       slopeStruct(ii).(sprintf('pId_%d',pId(jj))).ThetaFwd = tempData(jj,2);
       slopeStruct(ii).(sprintf('pId_%d',pId(jj))).ThetaRev = tempData(jj,3);
       slopeStruct(ii).(sprintf('pId_%d',pId(jj))).PhiFwd = tempData(jj,4);
       slopeStruct(ii).(sprintf('pId_%d',pId(jj))).PhiRev = tempData(jj,5);
    end
end

% generate plots in 3x3 grid, one for each positioner.
for jj=1:length(pId)
    
    pidStr = sprintf('pId_%d',pId(jj));
    % pull out data for one positioner
    td = [slopeStruct.(pidStr)];
    
    plotMe = [[td.ThetaFwd]; [td.ThetaRev]; [td.PhiFwd]; [td.PhiRev]].';
    
    plotMeDev  = bsxfun(@minus, plotMe, mean(plotMe));
    plotMeDev1 = bsxfun(@minus, plotMe, plotMe(1,:));

    %% slope vs. run
    figure(fh1)
    subplot(3,3,jj)
    plot(plotMe,'o-')
    title(pidStr,'Interpreter','none')
    xlabel('Target Run')
    ylabel('Err vs Request Slope')
    if jj == 1
        legend('ThetaFwd','ThetaRev','PhiFwd','PhiRev')
    end
        
    %% slope - <slope> vs. run
    figure(fh2)
    subplot(3,3,jj)
    plot(plotMeDev,'x-')
    title(pidStr,'Interpreter','none')
    xlabel('Target Run')
    ylabel('Deviation from Slope Mean')
    if jj == 1
        legend('ThetaFwd','ThetaRev','PhiFwd','PhiRev')
    end
    ylim([-.4,.4])

    %% slope - [first slope] vs. run

    figure(fh3)
    subplot(3,3,jj)
    plot(plotMeDev1,'x-')
    title(pidStr,'Interpreter','none')
    xlabel('Target Run')
    ylabel('Deviation from Initial Slope')
    if jj == 1
        legend('ThetaFwd','ThetaRev','PhiFwd','PhiRev')
    end
    ylim([-.4,.4])
end
