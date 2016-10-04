function output=analyzeTargetRun(pause_time)

warning('I changed this from a script to a function so that it wouldn''t pollute my base namespace -- Peter')

%% TARGET CONVERGENCE ANALYSIS SCRIPT
%clear all 
close all

%% INPUTS
% Provide logdirectory
logDir = '.';
saveFigs = false;
CCT = 5; %microns

if ~exist('pause_time','var'), pause_time = 0; end;
 

%% EXECUTION CODE (DANGER IF YOU MESS WITH THIS)

% $$$ [cfgFile,cfgPath] = uigetfile('*.xml','Select the original XML config file used for this test');
% $$$ CobraConfig = xml2struct(fullfile(cfgPath,cfgFile));
CobraConfig = loadCfgXml;

% Find the positioner target logs in current directory and store pids
posFiles = dir2cell('mId_*pId*.txt');
tempvar = regexp(posFiles,'(?<=mId_)\d*','match');
mId = cellfun(@str2num,[tempvar{:}]);
tempvar = regexp(posFiles,'(?<=pId_)\d*','match');
mpId = cellfun(@str2num,[tempvar{:}]);
clear tempvar;
pId = (mId-1)*57 + mpId;

% Initialize vars
cConvAll = [];
legendNames = {};

% For each log file
for ii=1:length(pId)
    pIdSTR = strcat('pId',num2str(pId(ii)));
    logfile = posFiles{ii};
    targetLogFile = fullfile(logDir,logfile);
    % Check target logfile for #INF instances and correct them if necessary
    checkLogFile(targetLogFile);
    
    dataName = char(regexp(logfile,'.*(?=\.txt)','match'));
    figPrefix = strcat(num2str(CCT),'um_',dataName);
    [tempData strArr] = processTargetLog(targetLogFile,pIdSTR,CobraConfig,CCT,saveFigs,figPrefix);
    cConvAll = [cConvAll; tempData.convP];
    legendNames = vertcat(legendNames,{dataName});
    eval([dataName ' = tempData;']);
    better_dataName = [regexprep(dataName,'(\d+)', ...
                                 '${sprintf(''%02d'',str2num($1))}') '_str'];
    eval([better_dataName ' = strArr;']);
    save(fullfile(logDir,strcat(num2str(CCT),'um_',dataName,'.mat')),dataName,better_dataName);
    drawnow;
    pause(pause_time);
    %s close all
   
end

% list the outputs you want here:
%output=packstruct(logDir,saveFigs);