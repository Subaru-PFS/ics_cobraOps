clear all;

%% Define test directory with folders containing the streak images
wsd = 'C:\Users\cmorantz\Dropbox\PFS_EM\TEST_RESULTS\StreakResults\DailyBenchmarkResults/';
%wsd = './';
file = 'workspace_streaks1';
wss = dir2cell([wsd '*.mat']);

dellist = [];
for jj = 1: length(wss)
    if(cellfun('length', strfind(wss(jj), file)) < 1)
        dellist = horzcat(dellist, jj);
    end
end
wss(dellist) = [];

for ps = 1:5
    pid = sprintf('pId%d',ps);
    meansFw.(pid) = [];
      meansRv.(pid) = [];
end


for jj = 1:length(wss)
    tmp = wss(jj);
    load(tmp{1},'meansityFw','meansityRv','motorMapFw','motorMapRv');
    for ps = 1:5
        pidfw = sprintf('pId%dfw',ps);
        pidrv = sprintf('pId%drv',ps);
        pidfwrsd = sprintf('pId%dfwrsd',ps);
        pidrvrsd = sprintf('pId%drvrsd',ps);
        pid = sprintf('pId%d',ps);
%         strkData(jj).(pidfw) = meansityFw.(pid);
%         strkData(jj).(pidrv) = meansityRv.(pid);
        strkData(jj).(pidfw) = motorMapFw.(pid)(2,:)';
        strkData(jj).(pidrv) = motorMapRv.(pid)(2,:)';
        strkData(jj).(pidfwrsd) = motorMapFw.(pid)(2,:)' / median(motorMapFw.(pid)(2,:));
        strkData(jj).(pidrvrsd) = motorMapRv.(pid)(2,:)' / median(motorMapRv.(pid)(2,:));
        meansFw.(pid) = vertcat(meansFw.(pid), meansityFw.(pid));
        meansRv.(pid) = vertcat(meansRv.(pid), meansityRv.(pid));
        
    end 
end


for ps = 1:5
    
    crange = [-1 1]*log10(2);
    
    fwfld = sprintf('pId%dfwrsd',ps);
    rvfld = sprintf('pId%drvrsd',ps);
    
    figure(100); set(gcf,'Name','Fwd strkData, log scale, normalized')
    subplot(2,3,ps)
    imagesc(log10([strkData.(fwfld)].'));
    caxis(crange);
    title(fwfld)
    xlabel('motor map bin')
    ylabel('motor map trial')
    
    figure(200); set(gcf,'Name','Rev strkData, log scale, normalized')
    subplot(2,3,ps)
    imagesc(log10([strkData.(rvfld)].'));
    caxis(crange);
    title(rvfld)
    xlabel('motor map bin')
    ylabel('motor map trial')
    
    crange = [0 .5];
    
    fwfld = sprintf('pId%dfw',ps);
    rvfld = sprintf('pId%drv',ps);
    
    figure(300); set(gcf,'Name','Fwd strkData, raw')
    subplot(2,3,ps)
    imagesc([strkData.(fwfld)].');
    caxis(crange);
    title(fwfld)
    xlabel('motor map bin')
    ylabel('motor map trial')
   
    figure(400); set(gcf,'Name','Rev strkData, raw')
    subplot(2,3,ps)
    imagesc([strkData.(rvfld)].');
    caxis(crange);
    title(rvfld)
    xlabel('motor map bin')
    ylabel('motor map trial')
end

% 
% 
% figure(1)
% imagesc(meansRv.pId1)
% %xlim([45 75])
% zlim([0 7*10^5])
% caxis([0 7*10^5])
% colorbar
% 
% 
% 
%  