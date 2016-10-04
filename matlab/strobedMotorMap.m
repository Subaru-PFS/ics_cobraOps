%% CLEAR WORKSPACE [OPTIONAL]
% clear all
% close all
% clc

%% INPUTS (CHANGE THIS SECTION TO SETUP MOTOR MAP)
% What is the title of this motormap. This will be the name of matlab structure containing motor map
title = 'Map1';

% What do the image files start with (don't change this for now)
FilterString = 'motorms';

% Name the image with the motor map strobe
mMapImage1 = 'motorms1fwi1_ImageId_1029';
mMapImage2 = 'motorms1fwi2_ImageId_1030';
% mMapImage1 = 'motorms2INref_ImageId_1033';
% mMapImage2 = 'motorms2fw_ImageId_1036';

title = strcat('thetaFwdMap_',char(regexp(mMapImage2,'ImageId_\d*','match')));

% Setup color for plots
CW = copper(3);
CWn = 3;

% saveDir = 'C:/Users/cmorantz/Dropbox/PFSMatlab/Results/GEN2_GEN3_Center_Metrology/'; %Use / and end with /

% Provide the directory containing the motor map image matlab files
% dataDir = ''; %Use \ 
% loadPrevious;
% Alternatively, load a specific mat file
datapath = fullfile(getenv('HOME'),'Dropbox','PFS_EM','TEST_RESULTS','Metrology','FiveEM_Archive');
load(fullfile(datapath,'strobedMaps_05_08_14_11_08_42.mat'));

% Load config file
% OLD: load('C:\Users\cmorantz\Dropbox\PFS_EM\SVN\MATLAB\EMconfig_03-31-14.mat');
% User selects XML config file
CobraConfig = loadCfgXml;

for pID = 1:9
  center(pID) = getARMval(CobraConfig,pID,'Global_base_pos_x') + ...
                getARMval(CobraConfig,pID,'Global_base_pos_y') * i;
  L1(pID)     = getARMval(CobraConfig,pID,'Link1_Link_Length');
  L2(pID)     = getARMval(CobraConfig,pID,'Link2_Link_Length');
end
%%
%%confirmed working to this point PHM 6/19/14
%%

% What direction is this motor map?
direction = 'FWD'; % Either REV or FWD

% Which stage is this for?
stage = 1;

% Provide the number of steps used for each move [GEN2,GEN3]
nstep = [200,200,200,200,200]; % S1 FWD 
% nstep = [150,150,150,150,150]; % S2 FWD 



%% GENERATE MOTOR MAP (DONT CHANGE THIS CODE BELOW)

% Evaluate first image
data1 = eval(mMapImage1);
data1.origin = 0;
% data1 = cobraCircle(data1,true);
assignin('base',mMapImage1,data1);
pids = regexp([fields(data1)],'pId\d+','match');
pids = pids(~cellfun('isempty',pids));

% if exist('mMapImage2','var')
    % Evaluate second image
    data2 = eval(mMapImage2);
    data2.origin = 0;
    assignin('base',mMapImage2,data2);

    % Create third combined image for full circle
    data3 = data1;
    for ii=1:length(pids)
        % Data 3 PID CCD pos' are combination of centroids from both images
        % but there is an additional row with either 1 or 2 value to
        % indicate which image. 
        data3.(char(pids{ii})).CCDpos = [data1.(char(pids{ii})).CCDpos, data2.(char(pids{ii})).CCDpos; 
                                            1*ones(1,length(data1.(char(pids{ii})).CCDpos)), 2*ones(1,length(data2.(char(pids{ii})).CCDpos))];
    end
% else
%     data3 = data1;
%     for ii=1:length(pids)
%         data3.(char(pids{ii})).CCDpos = [data1.(char(pids{ii})).CCDpos; 
%                                             1*ones(1,length(data1.(char(pids{ii})).CCDpos))];
%     end
% end

data3 = cobraCircle(data3,true);

tempMMap = motorMap(data3,nstep,stage,CobraConfig,direction);

assignin('base',title,tempMMap);




%% PLOT MAPS
fh = existingFigure('fwd motor maps');
for ii=1:5
    fldID = sprintf('pId%d',ii);
    subplot(2,3,ii)
    ph1 = plot(tempMMap.(fldID).mMap(1,:),tempMMap.(fldID).mMap(2,:),'r+-','color',CW(CWn,:));
    hl = appendLegend(ph1,strcat('05_08_14_11_08_42 StrobedMap',num2str(CWn)));
    set(hl,'Interpreter','none')
    hold on
%     axis auto
end
    

































