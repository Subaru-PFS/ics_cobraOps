%% TARGET CONVERGENCE ANALYSIS SCRIPT
clear all 
close all

%% INPUTS
% Provide arrays of corresponding module and positioner IDs. The same index
% in each array should match one positioner.
mId = [1,1,1,1,1,1,1];
pId = [1,2,3,4,7,8,9];
% Provide logdirectory
logDir = '..\..\TEST_RESULTS\TargetConvergence\06_17_14_11_21_32_TargetRun';
CobraConfigFile = 'CobraConfig_061014_5_inv.mat';
saveFigs = true;
CCT = 5; %microns



%% EXECUTION CODE (DANGER IF YOU MESS WITH THIS)
cConvAll = [];
legendNames = {};

for ii=1:length(pId)
    pIdSTR = strcat('pId',num2str(pId(ii)));
    logfile = strcat('mId_',num2str(mId(ii)),'_pId_',num2str(pId(ii)),'.txt');
    targetLogFile = fullfile(logDir,logfile);
    
    % Check target logfile for #INF instances and correct them if necessary
    checkLogFile(targetLogFile);
    
    dataName = char(regexp(logfile,'.*(?=\.txt)','match'));
    figPrefix = strcat(num2str(CCT),'um_',dataName);
    [tempData strArr] = processTargetLog(targetLogFile,pIdSTR,CobraConfigFile,CCT,saveFigs,figPrefix);
    cConvAll = [cConvAll; tempData.convP];
    legendNames = vertcat(legendNames,{dataName});
    assignin('base',dataName,tempData)
    assignin('base',strcat(dataName,'_str'),strArr)
    save(fullfile(logDir,strcat(num2str(CCT),'um_',dataName,'.mat')),dataName,strcat(dataName,'_str'))
%     pause;
    close all
end
 

