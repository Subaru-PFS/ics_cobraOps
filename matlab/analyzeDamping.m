clear all;
close all;

%% INPUTS
% Provide data directory
dataDir = '..\..\TEST_RESULTS\TargetConvergence\06_16_14_23_18_51_TargetRun';

% Provide cobra config structure file
cnfgFile = 'CobraConfig_061614_1.mat';

% Max steps
maxSteps = 16;

% Update motor maps?
revUpdate = false;
fwdUpdate = true;

% Is there  a motorMap workspace to compare with
% thetaMMapWorkspace = '..\..\TEST_RESULTS\StreakResults\workspacestage1.mat';
% phiMMapWorkspace = '..\..\TEST_RESULTS\StreakResults\workspacestage2.mat';

%% EXECUTION CODE
% Load all mat files in the dataDir
loadPrevious;

% Load cobra config structure
load(cnfgFile)

if exist('thetaMMapWorkspace','var'), 
    load(thetaMMapWorkspace)
    assignin('base','thetaMotorMapFw',motorMapFw)
    assignin('base','thetaMotorMapRv',motorMapRv)
end;

if exist('phiMMapWorkspace','var'), 
    load(phiMMapWorkspace)
    assignin('base','phiMotorMapFw',motorMapFw)
    assignin('base','phiMotorMapRv',motorMapRv)
end;

S = whos;
Names = {S.name};
SInd = find(~cellfun(@isempty,regexp(Names,'^mId_\d*_pId_\d*_str','match')));
c = 1;
fh1 = figure('Name','Error vs Request');
fh2 = figure('Name','Error Residuals');
CW = lines(length(SInd));
for ii=SInd
    pId = regexp(Names(ii),'(?<=pId_)\d*','match');
    pId = pId{1}{1};
    legendArr{c} = pId;
    cnfgID = ['ARM_DATA_' pId];
    
    qtemp = analyzeTC_CIT(eval(Names{ii}), CobraConfig, maxSteps);
    
    J1posSlope = qtemp.J1.posslope;
    J2posSlope = qtemp.J2.posslope;
    J1negSlope = qtemp.J1.negslope;
    J2negSlope = qtemp.J2.negslope;
    
    updateMaps = revUpdate || fwdUpdate;
    
    if updateMaps
        % Make a copy of original config
        CobraConfigOld = CobraConfig;
        
        % Get motor maps from config structure and convert to numeric
        % arrays
        
        Joint1_rev_stepsizes = mMapStr2NumArr(CobraConfig.ARM_DATA.(cnfgID).FAST_CALIBRATION_TABLE.Joint1_rev_stepsizes.Text);
        Joint2_rev_stepsizes = mMapStr2NumArr(CobraConfig.ARM_DATA.(cnfgID).FAST_CALIBRATION_TABLE.Joint2_rev_stepsizes.Text);
        Joint1_fwd_stepsizes = mMapStr2NumArr(CobraConfig.ARM_DATA.(cnfgID).FAST_CALIBRATION_TABLE.Joint1_fwd_stepsizes.Text);
        Joint2_fwd_stepsizes = mMapStr2NumArr(CobraConfig.ARM_DATA.(cnfgID).FAST_CALIBRATION_TABLE.Joint2_fwd_stepsizes.Text);
        

        % Initialize new motor maps as copies of original
        newJoint1_rev_stepsizes = Joint1_rev_stepsizes;
        newJoint2_rev_stepsizes = Joint2_rev_stepsizes;
        newJoint1_fwd_stepsizes = Joint1_fwd_stepsizes;
        newJoint2_fwd_stepsizes = Joint2_fwd_stepsizes;
        
        if fwdUpdate
            % Apply damping correction factor to motor map. First 2 entries are
            % size and number of entries in motor map.
            newJoint1_fwd_stepsizes(3:end) = Joint1_fwd_stepsizes(3:end) * (1-J1posSlope.a);
            newJoint2_fwd_stepsizes(3:end) = Joint2_fwd_stepsizes(3:end) * (1-J2posSlope.a);

            % Convert numeric motor map to comma-separated string and save to
            % structure
            CobraConfig.ARM_DATA.(cnfgID).FAST_CALIBRATION_TABLE.Joint1_fwd_stepsizes.Text = NumArr2mMapStr(newJoint1_fwd_stepsizes);
            CobraConfig.ARM_DATA.(cnfgID).FAST_CALIBRATION_TABLE.Joint2_fwd_stepsizes.Text = NumArr2mMapStr(newJoint2_fwd_stepsizes);

            % Copy FAST maps to SLOW maps
            CobraConfig.ARM_DATA.(cnfgID).SLOW_CALIBRATION_TABLE.Joint1_fwd_stepsizes.Text = CobraConfig.ARM_DATA.(cnfgID).FAST_CALIBRATION_TABLE.Joint1_fwd_stepsizes.Text;
            CobraConfig.ARM_DATA.(cnfgID).SLOW_CALIBRATION_TABLE.Joint2_fwd_stepsizes.Text = CobraConfig.ARM_DATA.(cnfgID).FAST_CALIBRATION_TABLE.Joint2_fwd_stepsizes.Text;
        end
        
        if revUpdate
            % Apply damping correction factor to motor map. First 2 entries are
            % size and number of entries in motor map.
            newJoint1_rev_stepsizes(3:end) = Joint1_rev_stepsizes(3:end) * (1-J1negSlope.a);
            newJoint2_rev_stepsizes(3:end) = Joint2_rev_stepsizes(3:end) * (1-J2negSlope.a);

            % Convert numeric motor map to comma-separated string and save to
            % structure
            CobraConfig.ARM_DATA.(cnfgID).FAST_CALIBRATION_TABLE.Joint1_rev_stepsizes.Text = NumArr2mMapStr(newJoint1_rev_stepsizes);
            CobraConfig.ARM_DATA.(cnfgID).FAST_CALIBRATION_TABLE.Joint2_rev_stepsizes.Text = NumArr2mMapStr(newJoint2_rev_stepsizes);

            % Copy FAST maps to SLOW maps
            CobraConfig.ARM_DATA.(cnfgID).SLOW_CALIBRATION_TABLE.Joint1_rev_stepsizes.Text = CobraConfig.ARM_DATA.(cnfgID).FAST_CALIBRATION_TABLE.Joint1_rev_stepsizes.Text;
            CobraConfig.ARM_DATA.(cnfgID).SLOW_CALIBRATION_TABLE.Joint2_rev_stepsizes.Text = CobraConfig.ARM_DATA.(cnfgID).FAST_CALIBRATION_TABLE.Joint2_rev_stepsizes.Text;
        end

    end
    
    figure(fh1)
    sph1 = subplot(211);
    ph1(c) = plot(qtemp.J1.req,qtemp.J1.err0,'.','color',CW(c,:));
    hold on
    plot(qtemp.J1.req,qtemp.J1.req*J1posSlope.a,'--','color',CW(c,:))
    sph2 = subplot(212);
    plot(qtemp.J2.req,qtemp.J2.err0,'.','color',CW(c,:))
    hold on
    plot(qtemp.J2.req,qtemp.J2.req*J2posSlope.a,'--','color',CW(c,:))
    
    figure(fh2)
    sph1 = subplot(211);
    plot(qtemp.J1.req,qtemp.J1.err,'.','color',CW(c,:))
    hold on
    sph2 = subplot(212);
    plot(qtemp.J2.req,qtemp.J2.err,'.','color',CW(c,:))
    hold on
    
    fhs(c) = figure('name',cnfgID);
    subplot(2,2,1)
    plot(newJoint1_fwd_stepsizes(3:end),'r')
    hold on
    plot(Joint1_fwd_stepsizes(3:end),'k')
    title('theta fwd')
    subplot(2,2,2)
    plot(newJoint1_rev_stepsizes(3:end),'r')
    hold on
    plot(Joint1_rev_stepsizes(3:end),'k')
    title('theta rev')
    subplot(2,2,3)
    plot(newJoint2_fwd_stepsizes(3:end),'r')
    hold on
    plot(Joint2_fwd_stepsizes(3:end),'k')
    title('phi fwd')
    subplot(2,2,4)
    plot(newJoint2_rev_stepsizes(3:end),'r')
    hold on
    plot(Joint2_rev_stepsizes(3:end),'k')
    title('phi rev')
    
    if exist('thetaMotorMapFw','var')
        pIdnum = str2num(pId);
        if pIdnum > 4
            pIdnum = pIdnum - 2;
        end
        pidStr = ['pId' num2str(pIdnum)];
        
        subplot(2,2,1)
        plot(thetaMotorMapFw.(pidStr)(2,3:end),'b--')
        subplot(2,2,2)
        plot(thetaMotorMapRv.(pidStr)(2,3:end),'b--')
        
        if length(thetaMotorMapFw.(pidStr)(2,:)) ~= 100
            keyboard
        end
        
        CobraConfigFromStreakMaps = CobraConfig;
        
        CobraConfigFromStreakMaps.ARM_DATA.(cnfgID).FAST_CALIBRATION_TABLE.Joint1_fwd_stepsizes.Text = ...
            NumArr2mMapStr([100,100,thetaMotorMapFw.(pidStr)(2,:)]);
        CobraConfigFromStreakMaps.ARM_DATA.(cnfgID).FAST_CALIBRATION_TABLE.Joint1_rev_stepsizes.Text = ...
            NumArr2mMapStr([100,100,thetaMotorMapRv.(pidStr)(2,:)]);
        CobraConfigFromStreakMaps.ARM_DATA.(cnfgID).SLOW_CALIBRATION_TABLE.Joint1_fwd_stepsizes.Text = ...
            CobraConfigFromStreakMaps.ARM_DATA.(cnfgID).FAST_CALIBRATION_TABLE.Joint1_fwd_stepsizes.Text;
        CobraConfigFromStreakMaps.ARM_DATA.(cnfgID).SLOW_CALIBRATION_TABLE.Joint1_rev_stepsizes.Text = ...
            CobraConfigFromStreakMaps.ARM_DATA.(cnfgID).FAST_CALIBRATION_TABLE.Joint1_rev_stepsizes.Text;
        
        CobraConfigFromStreakMaps.ARM_DATA.(cnfgID).FAST_CALIBRATION_TABLE.Joint1_fwd_regions.Text = ...
            NumArr2mMapStr([100,100,thetaMotorMapFw.(pidStr)(1,:)]);
        CobraConfigFromStreakMaps.ARM_DATA.(cnfgID).FAST_CALIBRATION_TABLE.Joint1_rev_regions.Text = ...
            NumArr2mMapStr([100,100,thetaMotorMapRv.(pidStr)(1,:)]);
        CobraConfigFromStreakMaps.ARM_DATA.(cnfgID).SLOW_CALIBRATION_TABLE.Joint1_fwd_regions.Text = ...
            CobraConfigFromStreakMaps.ARM_DATA.(cnfgID).FAST_CALIBRATION_TABLE.Joint1_fwd_regions.Text;
        CobraConfigFromStreakMaps.ARM_DATA.(cnfgID).SLOW_CALIBRATION_TABLE.Joint1_rev_regions.Text = ...
            CobraConfigFromStreakMaps.ARM_DATA.(cnfgID).FAST_CALIBRATION_TABLE.Joint1_rev_regions.Text;
    end
        
    if exist('phiMotorMapFw','var')
        pIdnum = str2num(pId);
        if pIdnum > 4
            pIdnum = pIdnum - 2;
        end
        pidStr = ['pId' num2str(pIdnum)];
        
        subplot(2,2,3)
        plot(phiMotorMapFw.(pidStr)(2,3:end),'b--')
        subplot(2,2,4)
        plot(phiMotorMapRv.(pidStr)(2,3:end),'b--')
        
        if length(phiMotorMapFw.(pidStr)(2,:)) ~= 45
            keyboard
        end
        
        CobraConfigFromStreakMaps = CobraConfig;
        
        CobraConfigFromStreakMaps.ARM_DATA.(cnfgID).FAST_CALIBRATION_TABLE.Joint2_fwd_stepsizes.Text = ...
            NumArr2mMapStr([45,100,phiMotorMapFw.(pidStr)(2,:)]);
        CobraConfigFromStreakMaps.ARM_DATA.(cnfgID).FAST_CALIBRATION_TABLE.Joint2_rev_stepsizes.Text = ...
            NumArr2mMapStr([45,100,phiMotorMapRv.(pidStr)(2,:)]);
        CobraConfigFromStreakMaps.ARM_DATA.(cnfgID).SLOW_CALIBRATION_TABLE.Joint2_fwd_stepsizes.Text = ...
            CobraConfigFromStreakMaps.ARM_DATA.(cnfgID).FAST_CALIBRATION_TABLE.Joint2_fwd_stepsizes.Text;
        CobraConfigFromStreakMaps.ARM_DATA.(cnfgID).SLOW_CALIBRATION_TABLE.Joint2_rev_stepsizes.Text = ...
            CobraConfigFromStreakMaps.ARM_DATA.(cnfgID).FAST_CALIBRATION_TABLE.Joint2_rev_stepsizes.Text;
        
        CobraConfigFromStreakMaps.ARM_DATA.(cnfgID).FAST_CALIBRATION_TABLE.Joint2_fwd_regions.Text = ...
            NumArr2mMapStr([45,100,phiMotorMapFw.(pidStr)(1,:)]);
        CobraConfigFromStreakMaps.ARM_DATA.(cnfgID).FAST_CALIBRATION_TABLE.Joint2_rev_regions.Text = ...
            NumArr2mMapStr([45,100,phiMotorMapRv.(pidStr)(1,:)]);
        CobraConfigFromStreakMaps.ARM_DATA.(cnfgID).SLOW_CALIBRATION_TABLE.Joint2_fwd_regions.Text = ...
            CobraConfigFromStreakMaps.ARM_DATA.(cnfgID).FAST_CALIBRATION_TABLE.Joint2_fwd_regions.Text;
        CobraConfigFromStreakMaps.ARM_DATA.(cnfgID).SLOW_CALIBRATION_TABLE.Joint2_rev_regions.Text = ...
            CobraConfigFromStreakMaps.ARM_DATA.(cnfgID).FAST_CALIBRATION_TABLE.Joint2_rev_regions.Text;
    end
    
    assignin('base',sprintf('q%d',c),qtemp);
    c=c+1;
end

% Assign cobra configs more reasonable names
CobraConfigAdj = CobraConfig;
clear CobraConfig

figure(fh1)

subplot(211)
title('\Theta');    
xlabel('request angle [rad]');
ylabel('error angle [rad]');
refline(0,0,0,'k');
refline(0,0,Inf,'k');
text(0,1,'  OVERDAMPED');
text(0,-.5,'  UNDERDAMPED');
text(0,1,'UNDERDAMPED  ','horiz','right');
text(0,-.5,'OVERDAMPED  ','horiz','right');

subplot(212)
title('\Phi');
xlabel('request angle [rad]');
ylabel('error angle [rad]');
refline(0,0,0,'k');
refline(0,0,Inf,'k');
text(0,.1,'  OVERDAMPED');
text(0,-.1,'  UNDERDAMPED');
text(0,.1,'UNDERDAMPED  ','horiz','right');
text(0,-.1,'OVERDAMPED  ','horiz','right');

legend(ph1,[legendArr],'location','NW')

set(findall(gcf,'type','text'),'fontSize',14)

saveas(gcf,'error_vs_request_plots','fig');

figure(fh2);
subplot(211)
hold off;
title('\Theta - corrected');
xlabel('request angle [rad]');
ylabel('error angle [rad]');
refline(0,0,0,'k');
refline(0,0,Inf,'k');
text(0,1,'  OVERDAMPED');
text(0,-.5,'  UNDERDAMPED');
text(0,1,'UNDERDAMPED  ','horiz','right');
text(0,-.5,'OVERDAMPED  ','horiz','right');
subplot(212)
hold off;
title('\Phi - corrected');
xlabel('request angle [rad]');
ylabel('error angle [rad]');
refline(0,0,0,'k');
refline(0,0,Inf,'k');
text(0,.1,'  OVERDAMPED');
text(0,-.1,'  UNDERDAMPED');
text(0,.1,'UNDERDAMPED  ','horiz','right');
text(0,-.1,'OVERDAMPED  ','horiz','right');
legend([legendArr],'location','NW')
set(findall(gcf,'type','text'),'fontSize',14)
drawnow;        
        










    
