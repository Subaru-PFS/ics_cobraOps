%% Convergence Testing Comparison Script
% This script takes the data directories you provide it which have mat
% files from previously run targConvAnlys script and compares the
% cummulative convergences of the different tests which those directories
% correspond to.
clear all
close all

%% INPUTS
dataDirs = {'06_17_14_22_41_35_TargetRun';
      '07_09_14_16_37_46_TargetRun'};
%     '06_18_14_21_43_16_TargetRun';
%     '06_19_14_13_44_30_TargetRun'};

loadFilter = '5um';

%% 
nTests = length(dataDirs);
for ii=1:nTests
    dataDir = dataDirs{1};
    loadPrevious
    S = whos;
    Names = {S.name};
    posFileInd = strmatch('mId',Names)';
    for vv = posFileInd
        vname = char(Names(vv));
        assignin('base',strcat('D',num2str(ii),'_',vname),eval(vname));
        clear (vname);
    end
end
    


% c=c+1;
% dataDir = 'C:\Users\cmorantz\Dropbox\PFS_EM\TEST_RESULTS\TargetConvergence\04_02_14_11_47_52_Targets_NoDamping';
% dscr{c} = strcat('D',num2str(c),' = Map1');
% loadPrevious
% S = whos;
% Names = {S.name};
% for vv = strmatch('mId',Names)'
%     vname = char(Names(vv));
%     assignin('base',strcat('D',num2str(c),'_',vname),eval(vname));
%     clear (vname);
% end
% % 
% c=c+1;
% dataDir = 'C:\Users\cmorantz\Dropbox\PFS_EM\TEST_RESULTS\TargetConvergence\04_02_14_17_49_31_CustomDamping';
% dscr{c} = strcat('D',num2str(c),' = Map1 w/ Estimated Adjustments');
% loadPrevious
% S = whos;
% Names = {S.name};
% for vv = strmatch('mId',Names)'
%     vname = char(Names(vv));
%     assignin('base',strcat('D',num2str(c),'_',vname),eval(vname));
%     clear (vname);
% end
% % 
% c=c+1;
% dataDir = 'C:\Users\cmorantz\Dropbox\PFS_EM\TEST_RESULTS\TargetConvergence\04_03_14_16_15_49_Targets';
% dscr{c} = strcat('D',num2str(c),' = Map1 w/ Calculated Adjustments');
% loadPrevious
% S = whos;
% Names = {S.name};
% for vv = strmatch('mId',Names)'
%     vname = char(Names(vv));
%     assignin('base',strcat('D',num2str(c),'_',vname),eval(vname));
%     clear (vname);
% end
% % 
% c=c+1;
% dataDir = 'C:\Users\cmorantz\Dropbox\PFS_EM\TEST_RESULTS\TargetConvergence\04_04_14_20_23_02_Targets';
% dscr{c} = strcat('D',num2str(c),' = Map1 Fixed PID3');
% loadPrevious
% S = whos;
% Names = {S.name};
% for vv = strmatch('mId',Names)'
%     vname = char(Names(vv));
%     assignin('base',strcat('D',num2str(c),'_',vname),eval(vname));
%     clear (vname);
% end
% % 
% c=c+1;
% dataDir = 'C:\Users\cmorantz\Dropbox\PFS_EM\TEST_RESULTS\TargetConvergence\04_05_14_18_36_53_Targets';
% dscr{c} = strcat('D',num2str(c),' = Map1 w/ More Adjustements');
% loadPrevious
% S = whos;
% Names = {S.name};
% for vv = strmatch('mId',Names)'
%     vname = char(Names(vv));
%     assignin('base',strcat('D',num2str(c),'_',vname),eval(vname));
%     clear (vname);
% end
% % 
% c=c+1;
% dataDir = 'C:\Users\cmorantz\Dropbox\PFS_EM\TEST_RESULTS\TargetConvergence\04_07_14_15_41_24_Targets';
% dscr{c} = strcat('D',num2str(c),' = Map1 further PID5 Adjustments');
% loadPrevious
% S = whos;
% Names = {S.name};
% for vv = strmatch('mId',Names)'
%     vname = char(Names(vv));
%     assignin('base',strcat('D',num2str(c),'_',vname),eval(vname));
%     clear (vname);
% end
% % 
% % c=c+1;
% % dataDir = 'C:\Users\cmorantz\Dropbox\PFS_EM\TEST_RESULTS\TargetConvergence\04_09_14_18_06_55_Targets';
% % dscr{c} = strcat('D',num2str(c),' = Adjusted map 04_09_14_18_06_55');
% % loadPrevious
% % S = whos;
% % Names = {S.name};
% % for vv = strmatch('mId',Names)'
% %     vname = char(Names(vv));
% %     assignin('base',strcat('D',num2str(c),'_',vname),eval(vname));
% %     clear (vname);
% % end
% % 
% c=c+1;
% dataDir = 'C:\Users\cmorantz\Dropbox\PFS_EM\TEST_RESULTS\TargetConvergence\04_10_14_14_30_50_Targets';
% dscr{c} = strcat('D',num2str(c),' = Map2');
% loadPrevious
% S = whos;
% Names = {S.name};
% for vv = strmatch('mId',Names)'
%     vname = char(Names(vv));
%     assignin('base',strcat('D',num2str(c),'_',vname),eval(vname));
%     clear (vname);
% end
% 
% c=c+1;
% dataDir = 'C:\Users\cmorantz\Dropbox\PFS_EM\TEST_RESULTS\TargetConvergence\04_13_14_07_42_20_Targets';
% dscr{c} = strcat('D',num2str(c),' = Map3');
% loadPrevious
% S = whos;
% Names = {S.name};
% for vv = strmatch('mId',Names)'
%     vname = char(Names(vv));
%     assignin('base',strcat('D',num2str(c),'_',vname),eval(vname));
%     clear (vname);
% end

% c=c+1;
% dataDir = 'C:\Users\cmorantz\Dropbox\PFS_EM\TEST_RESULTS\TargetConvergence\04_16_14_08_56_36_Targets';
% dscr{c} = strcat('D',num2str(c),' = Map4 for PID3 Phi');
% loadPrevious
% S = whos;
% Names = {S.name};
% for vv = strmatch('mId',Names)'
%     vname = char(Names(vv));
%     assignin('base',strcat('D',num2str(c),'_',vname),eval(vname));
%     clear (vname);
% end
% 
% 
% % 
% c=c+1;
% dataDir = 'C:\Users\cmorantz\Dropbox\PFS_EM\TEST_RESULTS\TargetConvergence\04_16_14_15_20_40_Targets';
% dscr{c} = strcat('D',num2str(c),' = PID3 fast/slow');
% loadPrevious
% S = whos;
% Names = {S.name};
% for vv = strmatch('mId',Names)'
%     vname = char(Names(vv));
%     assignin('base',strcat('D',num2str(c),'_',vname),eval(vname));
%     clear (vname);
% end
% 
% c=c+1;
% dataDir = 'C:\Users\cmorantz\Dropbox\PFS_EM\TEST_RESULTS\TargetConvergence\04_20_14_10_44_35_Targets';
% dscr{c} = strcat('D',num2str(c),' = PID3 fast/slow adjusted');
% loadPrevious
% S = whos;
% Names = {S.name};
% for vv = strmatch('mId',Names)'
%     vname = char(Names(vv));
%     assignin('base',strcat('D',num2str(c),'_',vname),eval(vname));
%     clear (vname);
% end
% 
% c=c+1;
% dataDir = 'C:\Users\cmorantz\Dropbox\PFS_EM\TEST_RESULTS\TargetConvergence\04_27_14_21_01_19_Targets';
% dscr{c} = strcat('D',num2str(c),' = Slow Map Updating');
% loadPrevious
% S = whos;
% Names = {S.name};
% for vv = strmatch('mId',Names)'
%     vname = char(Names(vv));
%     assignin('base',strcat('D',num2str(c),'_',vname),eval(vname));
%     clear (vname);
% end
% 
% c=c+1;
% dataDir = 'C:\Users\cmorantz\Dropbox\PFS_EM\TEST_RESULTS\TargetConvergence\04_28_14_11_13_33_Targets';
% dscr{c} = strcat('D',num2str(c),' = Test13: SlowUpdates WeightAvg & DampAdj');
% loadPrevious
% S = whos;
% Names = {S.name};
% for vv = strmatch('mId',Names)'
%     vname = char(Names(vv));
%     assignin('base',strcat('D',num2str(c),'_',vname),eval(vname));
%     clear (vname);
% end

% c=c+1;
% dataDir = 'C:\Users\cmorantz\Dropbox\PFS_EM\TEST_RESULTS\TargetConvergence\04_29_14_21_58_02_Targets';
% dscr{c} = strcat('D',num2str(c),' = Map4');
% loadPrevious
% S = whos;
% Names = {S.name};
% for vv = strmatch('mId',Names)'
%     vname = char(Names(vv));
%     assignin('base',strcat('D',num2str(c),'_',vname),eval(vname));
%     clear (vname);
% end

% c=c+1;
% dataDir = 'C:\Users\cmorantz\Dropbox\PFS_EM\TEST_RESULTS\TargetConvergence\DotObscureTargets\04_30_14_16_30_06';
% dscr{c} = strcat('D',num2str(c),' = Map4');
% loadPrevious
% S = whos;
% Names = {S.name};
% for vv = strmatch('mId',Names)'
%     vname = char(Names(vv));
%     assignin('base',strcat('D',num2str(c),'_',vname),eval(vname));
%     clear (vname);
% end

% c=c+1;
% dataDir = 'C:\Users\cmorantz\Dropbox\PFS_EM\TEST_RESULTS\TargetConvergence\05_06_14_14_36_34_Targets';
% dscr{c} = strcat('D',num2str(c),' = Test16: SlowUpdates StraightAvg & DampAdj');
% loadPrevious
% S = whos;
% Names = {S.name};
% for vv = strmatch('mId',Names)'
%     vname = char(Names(vv));
%     assignin('base',strcat('D',num2str(c),'_',vname),eval(vname));
%     clear (vname);
% end
% 
% c=c+1;
% dataDir = 'C:\Users\cmorantz\Dropbox\PFS_EM\TEST_RESULTS\TargetConvergence\05_07_14_20_25_12_Targets';
% dscr{c} = strcat('D',num2str(c),' = Test17: StreakMap5');
% loadPrevious
% S = whos;
% Names = {S.name};
% for vv = strmatch('mId',Names)'
%     vname = char(Names(vv));
%     assignin('base',strcat('D',num2str(c),'_',vname),eval(vname));
%     clear (vname);
% end

% c=c+1;
% dataDir = '..\..\TEST_RESULTS\TargetConvergence\06_12_14_17_42_30_TargetRun';
% dscr{c} = strcat('D',num2str(c),' = StrkMaps1_1');
% loadPrevious
% S = whos;
% Names = {S.name};
% posFileInd = strmatch('mId',Names)';
% for vv = posFileInd
%     vname = char(Names(vv))
%     assignin('base',strcat('D',num2str(c),'_',vname),eval(vname));
%     clear (vname);
% end
% 
% c=c+1;
% dataDir = '..\..\TEST_RESULTS\TargetConvergence\06_16_14_10_36_21_TargetRun';
% dscr{c} = strcat('D',num2str(c),' = StrkMaps1_2');
% loadPrevious
% S = whos;
% Names = {S.name};
% posFileInd = strmatch('mId',Names)';
% for vv = posFileInd
%     vname = char(Names(vv))
%     assignin('base',strcat('D',num2str(c),'_',vname),eval(vname));
%     clear (vname);
% end
% 
% c=c+1;
% dataDir = '..\..\TEST_RESULTS\TargetConvergence\06_16_14_23_18_51_TargetRun';
% dscr{c} = strcat('D',num2str(c),' = StrkMaps2_1');
% loadPrevious
% S = whos;
% Names = {S.name};
% posFileInd = strmatch('mId',Names)';
% for vv = posFileInd
%     vname = char(Names(vv))
%     assignin('base',strcat('D',num2str(c),'_',vname),eval(vname));
%     clear (vname);
% end



S = whos;
Names = {S.name};
SInd = find(~cellfun(@isempty,regexp(Names,'^D\d*_mId_\d*_pId_\d*$','match')));
cellpids = regexp([Names{:}],'(?<=D\d*_mId_\d*_pId_)\d*','match');
numpids = str2num(char(cellpids));
Npos = length(unique(numpids));
CW = colormap(hsv(nTests));
t=1;
ph = [];
% lgnd = [];
fh1 = figure;
conv7M = [];
nRows = floor(Npos/2);
nCol = ceil(Npos/nRows);
nIter = 0;
for ii = SInd
    figure(fh1);
    thisvar = eval(Names{ii});
    D = str2num(char(regexp(Names{ii},'(?<=D)\d*(?=_mId)','match')));
    id = str2num(char(regexp(Names{ii},'(?<=pId_)\d*','match')));
    sph = subplot(nRows,nCol,id);
    ph = plot(thisvar.convP,'Color',CW(D,:));
    tnIter = length(thisvar.convP);
    if tnIter > nIter, nIter=tnIter, end;
    lh1 = appendLegend(ph,Names{ii},'Interpreter','none');
    set(lh1,'Interpreter','none')
%     lgnd = [lgnd; Names{ii}];
    hold on
    t = t+1;
    title(strcat('PID',num2str(id),' 5um Cumulative Convergence'))
    ylabel('Percentage Convergence (%)')
    xlabel('Iterations')
    set(gca,'YTick',linspace(0,100,11))
    set(gca,'XTick',linspace(1,15,15))
    grid on
    ylim([0 100])
    l2 = plot(95*ones(1,nIter),'--g');
%     lh2 = appendLegend(l2,'95%');
%     set(lh2,'Interpreter','none')
    
    conv7M(D,id) = thisvar.convP(7);
    
    ccfh = figure('name',strcat('PID',num2str(id),' 5um Cumulative Convergence'));
%     ccfhs1 = subplot(1,1,1);
    ccfhs1 = axes;
    copyobj(allchild(sph),ccfhs1);
    copyobj(lh1,ccfh);
    ylabel('Percentage Convergence (%)')
    xlabel('Iterations')
    set(gca,'YTick',linspace(0,100,11))
    set(gca,'XTick',linspace(1,15,15))
    ylim([0 100])
    grid on
    
end


% Create textbox
th1 = annotation(fh1,'textbox',...
    [0.70 0.40 0.15 0.15],...
    'String',dscr);
set(th1, 'Interpreter','none')




% 
% S = whos;
% Names = {S.name};
% SInd = find(~cellfun(@isempty,regexp(Names,'^D\d*_mId_\d*_pId_\d*_str$','match')));
% CW = colormap(hsv(length(SInd)));
% t=1;
% ph = [];
% lgnd = [];
% fh2 = figure('Name','Distance to target');
% fh3 = figure('Name','J1 Errors');
% fh4 = figure('Name','J2 Errors')
% for ii = SInd
%     thisvar = eval(Names{ii});
%     id = str2num(char(regexp(Names{ii},'(?<=pId_)\d*','match')));
%     
%     figure(fh2)
%     subplot(nRows,nCol,id)
%     plot(log([thisvar.dist]));
%     hold on
%     
%     figure(fh3)
%     subplot(nRows,nCol,id)
%     plot([thisvar.J1err]);
%     xlabel('Iteration')
%     ylabel('J1 Error')
%     title(strcat('pid',num2str(id),' J1 Errors'))
%     
%     figure(fh4)
%     subplot(nRows,nCol,id)
%     plot([thisvar.J2err]);
%     xlabel('Iteration')
%     ylabel('J2 Error')
%     title(strcat('pid',num2str(id),' J2 Errors'))
%     
%     t = t+1;
% 
%     
%     
%     
%     
% end
% lh = legend(ph,lgnd);
% set(lh, 'Interpreter','none')