%Plot motor maps
%% Create Motor Map from HS
function output=plotMotMaps()


%% THETA MAPS
Files = {
%     'C:\Users\cmorantz\Dropbox\PFS_EM\TEST_RESULTS\StreakResults\Stages1\03_31_14\workspacestage1.mat';
%     'C:\Users\cmorantz\Dropbox\PFS_EM\TEST_RESULTS\StreakResults\Stages1\04_09_14_15_47_42\workspacestage1.mat';
    'C:\Users\cmorantz\Dropbox\PFS_EM\TEST_RESULTS\StreakResults\Stages1\04_11_14_11_18_21\workspacestage1.mat';
    'C:\Users\cmorantz\Dropbox\PFS_EM\TEST_RESULTS\StreakResults\Stages1\04_28_14_20_06_17_Streaks\workspacestage1_2.mat';
    'C:\Users\cmorantz\Dropbox\PFS_EM\TEST_RESULTS\StreakResults\Stages1\05_07_14_09_57_34\workspacestage1.mat';
    };
fwOnTimes = [
%     .07,.07,.05,.09,.2; % 3-30-14 restore point of 2014_03_27_StreaksForMM01.lst. Also confirmed in D:\PfsTests\Xml\bak-04-01-2014-EM_5Cobras target config
%     .07,.07,.05,.09,.3; % From 2014_04_09_StreaksForpId5test.lst stamped 4-9-14@15:47
    .09,.09,.07,.11,.3; % from 4/11/14 stamped version of 2014_03_27_StreaksForMM01.lst
    .09,.09,.07,.11,.3;  % 2014_04_28_StreaksForMM01 stamped 4-28-14@20:06
    .09,.09,.07,.11,.3; % 2014_05_07_StreaksForMMIncludingDots01.lst
    ];
    
rvOnTimes = [
%     .1,.11,.05,.16,.2;  % 3-30-14 restore point of 2014_03_27_StreaksForMM01.lst. Also confirmed in D:\PfsTests\Xml\bak-04-01-2014-EM_5Cobras target config
%     .1,.11,.05,.16,.3;  % From 2014_04_09_StreaksForpId5test.lst stamped 4-9-14@15:47
    .12,.13,.07,.18,.3; % from 4/11/14 stamped version of 2014_03_27_StreaksForMM01.lst
    .12,.13,.07,.18,.3;  % 2014_04_28_StreaksForMM01 stamped 4-28-14@20:06
    .12,.13,.07,.18,.3; % 2014_05_07_StreaksForMMIncludingDots01.lst
    ];


%% PHI MAPS
% Files = {
% %     'C:\Users\cmorantz\Dropbox\PFS_EM\TEST_RESULTS\StreakResults\Stages2\03_31_14\workspacestage2.mat';
% %     'C:\Users\cmorantz\Dropbox\PFS_EM\TEST_RESULTS\StreakResults\Stages2\04_11_14_11_18_21\workspacestage2.mat';
%     'C:\Users\cmorantz\Dropbox\PFS_EM\TEST_RESULTS\StreakResults\Stages2\04_15_14_14_56_13\workspacestage2.mat';
%     'C:\Users\cmorantz\Dropbox\PFS_EM\TEST_RESULTS\StreakResults\Stages2\04_28_14_20_06_17_Streaks\workspacestage2.mat';
%     'C:\Users\cmorantz\Dropbox\PFS_EM\TEST_RESULTS\StreakResults\Stages2\05_07_14_09_57_34\workspacestage2.mat';
%     };
% fwOnTimes = [
% %     .05,.04,.15,.04,.03;    %3-30-14 restore point of 2014_03_27_StreaksForMM01.lst
% %     .06,.05,.16,.05,.04;    % from 4/11/14 stamped version of 2014_03_27_StreaksForMM01.lst
%     .06,.05,.06,.05,.04;    % from 2014_04_15_StreaksForMM01.lst stamped 4-15-14 @ 14:55
%     .06,.05,.06,.05,.04;     % 2014_04_28_StreaksForMM01 stamped 4-28-14@20:06
%     .06,.05,.06,.05,.04;
%     ];
% 
% rvOnTimes = [
% %     .05,.06,.18,.04,.04;    %3-30-14 restore point of 2014_03_27_StreaksForMM01.lst. Also confirmed in D:\PfsTests\Xml\bak-04-01-2014-EM_5Cobras target config
% %     .06,.07,.19,.05,.05;    % from 4/11/14 stamped version of 2014_03_27_StreaksForMM01.lst
%     .06,.07,.07,.05,.05;    % from 2014_04_15_StreaksForMM01.lst stamped 4-15-14 @ 14:55
%     .06,.07,.07,.05,.05;     % 2014_04_28_StreaksForMM01 stamped 4-28-14@20:06
%     .06,.07,.07,.05,.05;
%     ];


%% EXECTUTION CODE

CW = colormap(lines);
fwdFig = figure('name','fwd motor maps');
rvFig = figure('name','rev motor maps');
for ii = 1:length(Files)
    name = char(regexp(Files{ii},'(?<=Stages\d\\)[\d_-]*','match'));
    plotmap(Files{ii},CW(ii,:),name,fwOnTimes(ii,:),rvOnTimes(ii,:));
end






end


%% Plot Map Function
function output = plotmap(mapfile,colorspec,name,fwOn,rvOn)
load(mapfile)
for ps = [1:length(crm(:,1))]
    pid = sprintf('pId%d',ps) ; 
    
    %Plot fwd maps
    fwdFig = existingFigure('fwd motor maps');
    subplot(2,3,ps)
    tp = plot(motorMapFw.(pid)(1,:),motorMapFw.(pid)(2,:), '-x','color',colorspec);
        xlabel('Joint angle')
    ylabel('deg/step')
    title(strcat(pid,' ',stage,' motor map'))
    tl = appendLegend(tp,strcat(name,' dur=',num2str(fwOn(ps)))); 
    set(tl,'Interpreter','none')
    if stage == 'stage2'
        xlim([0,200])
        ylim([0,1])
    elseif stage == 'stage1'
        xlim([0,400])
        ylim([0,1])
    end
    hold on
    
    
    
    % Plot rev maps
    rvFig = existingFigure('rev motor maps');
    subplot(2,3,ps)
    tp = plot(motorMapRv.(pid)(1,:),motorMapRv.(pid)(2,:), '--x','color',colorspec);
    xlabel('Joint angle')
    ylabel('deg/step')
    title(strcat(pid,' ',stage,' motor map'))
%     tl = appendLegend(tp,name); 
    tl = appendLegend(tp,strcat(name,' dur=',num2str(rvOn(ps))));
    set(tl,'Interpreter','none')
    if stage == 'stage2'
        xlim([0,200])
        ylim([0,1])
    elseif stage == 'stage1'
        xlim([0,400])
        ylim([0,1])
    end
    hold on
    
    output = 0;

end
end
