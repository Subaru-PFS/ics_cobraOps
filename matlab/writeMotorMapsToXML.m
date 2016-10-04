%Specifiy XML file
clear all
close all

% Choose xml file to start with
CobraConfig = loadCfgXml;

xmlEnd = fullfile(CobraConfig.cfgPath,['new_' CobraConfig.cfgFile]);

[thetaWorkspace thetapath] = uigetfile('D:\PfsTests\*.mat','Select theta streak map workspace');
load(fullfile(thetapath,thetaWorkspace))
motorMapJ1Fw = motorMapFw;
motorMapJ1Rv = motorMapRv;
[phiWorkspace phipath] = uigetfile('D:\PfsTests\*.mat','Select phi streak map workspace');
load(fullfile(phipath,phiWorkspace))
motorMapJ2Rv = motorMapRv;
motorMapJ2Fw = motorMapFw;

keyboard;

J1fwOnTimes = [.07,.07,.13,.07,.18,.25,.07];
J1rvOnTimes = [.13,.09,.13,.07,.18,.3,.09];
J2fwOnTimes = [.05,.05,.05,.05,.05,.05,.05];
J2rvOnTimes = [.05,.05,.05,.05,.05,.05,.05];

J1fwIntTimes = [2.5,2.5,2.5,2.5,2.5,2.5,2.5];
J1rvIntTimes = [2.5,2.5,2.5,2.5,2.5,2.5,2.5];
J2fwIntTimes = [2.5,2.5,2.5,2.5,2.5,2.5,2.5];
J2rvIntTimes = [2.5,2.5,2.5,2.5,2.5,2.5,2.5];

for pid=1:7
    
    if pid > 4
        cpid = pid+2;
    else
        cpid = pid;
    end
    
    fldID = sprintf('pId%d',pid);
    cnfgID = sprintf('ARM_DATA_%d',cpid);
    
    CobraConfig.ARM_DATA.(cnfgID).KINEMATICS.Link1_fwd_Duration = J1fwOnTimes(pid);
    CobraConfig.ARM_DATA.(cnfgID).KINEMATICS.Link1_fwd_Duration_Slow = J1fwOnTimes(pid); 
    CobraConfig.ARM_DATA.(cnfgID).KINEMATICS.Link1_fwd_Intervals = J1fwIntTimes(pid);
    CobraConfig.ARM_DATA.(cnfgID).KINEMATICS.Link1_fwd_Intervals_Slow = J1fwIntTimes(pid);
    CobraConfig.ARM_DATA.(cnfgID).KINEMATICS.Link1_rev_Duration = J1rvOnTimes(pid);
    CobraConfig.ARM_DATA.(cnfgID).KINEMATICS.Link1_rev_Duration_Slow = J1rvOnTimes(pid);
    CobraConfig.ARM_DATA.(cnfgID).KINEMATICS.Link1_rev_Intervals = J1rvIntTimes(pid);
    CobraConfig.ARM_DATA.(cnfgID).KINEMATICS.Link1_rev_Intervals_Slow = J1rvIntTimes(pid);
    
    CobraConfig.ARM_DATA.(cnfgID).KINEMATICS.Link2_fwd_Duration = J2fwOnTimes(pid);
    CobraConfig.ARM_DATA.(cnfgID).KINEMATICS.Link2_fwd_Duration_Slow = J2fwOnTimes(pid);
    CobraConfig.ARM_DATA.(cnfgID).KINEMATICS.Link2_fwd_Intervals = J2fwIntTimes(pid);
    CobraConfig.ARM_DATA.(cnfgID).KINEMATICS.Link2_fwd_Intervals_Slow = J2fwIntTimes(pid);
    CobraConfig.ARM_DATA.(cnfgID).KINEMATICS.Link2_rev_Duration = J2rvOnTimes(pid);
    CobraConfig.ARM_DATA.(cnfgID).KINEMATICS.Link2_rev_Duration_Slow = J2rvOnTimes(pid);
    CobraConfig.ARM_DATA.(cnfgID).KINEMATICS.Link2_rev_Intervals = J2rvIntTimes(pid);
    CobraConfig.ARM_DATA.(cnfgID).KINEMATICS.Link2_rev_Intervals_Slow = J2rvIntTimes(pid);
    
    CobraConfig.ARM_DATA.(cnfgID).KINEMATICS.Joint1_transition_angle = 400;
    CobraConfig.ARM_DATA.(cnfgID).KINEMATICS.Joint2_transition_angle = 400;
    
    %%FAST MOTOR MAPS
    %Joint 1 FW
    CobraConfig.ARM_DATA.(cnfgID).FAST_CALIBRATION_TABLE.Joint1_fwd_regions.Text =  array2text(motorMapJ1Fw.(fldID)(1,:));
    CobraConfig.ARM_DATA.(cnfgID).FAST_CALIBRATION_TABLE.Joint1_fwd_stepsizes.Text = array2text(motorMapJ1Fw.(fldID)(2,:));
    % Joint 2 FW
    CobraConfig.ARM_DATA.(cnfgID).FAST_CALIBRATION_TABLE.Joint2_fwd_regions.Text = array2text(motorMapJ2Fw.(fldID)(1,:));
    CobraConfig.ARM_DATA.(cnfgID).FAST_CALIBRATION_TABLE.Joint2_fwd_stepsizes.Text = array2text(motorMapJ2Fw.(fldID)(2,:));

    % Joint 1 RV
    CobraConfig.ARM_DATA.(cnfgID).FAST_CALIBRATION_TABLE.Joint1_rev_regions.Text =  array2text(motorMapJ1Rv.(fldID)(1,:));
    CobraConfig.ARM_DATA.(cnfgID).FAST_CALIBRATION_TABLE.Joint1_rev_stepsizes.Text = array2text(motorMapJ1Rv.(fldID)(2,:));
    % Joint 2 RV
    CobraConfig.ARM_DATA.(cnfgID).FAST_CALIBRATION_TABLE.Joint2_rev_regions.Text = array2text(motorMapJ2Rv.(fldID)(1,:));
    CobraConfig.ARM_DATA.(cnfgID).FAST_CALIBRATION_TABLE.Joint2_rev_stepsizes.Text = array2text(motorMapJ2Rv.(fldID)(2,:));


    %% SLOW MOTOR MAPS
    %Joint 1 FW
    CobraConfig.ARM_DATA.(cnfgID).SLOW_CALIBRATION_TABLE.Joint1_fwd_regions.Text =  array2text(motorMapJ1Fw.(fldID)(1,:));
    CobraConfig.ARM_DATA.(cnfgID).SLOW_CALIBRATION_TABLE.Joint1_fwd_stepsizes.Text = array2text(motorMapJ1Fw.(fldID)(2,:));
    % Joint 2 FW
    CobraConfig.ARM_DATA.(cnfgID).SLOW_CALIBRATION_TABLE.Joint2_fwd_regions.Text = array2text(motorMapJ2Fw.(fldID)(1,:));
    CobraConfig.ARM_DATA.(cnfgID).SLOW_CALIBRATION_TABLE.Joint2_fwd_stepsizes.Text = array2text(motorMapJ2Fw.(fldID)(2,:));

    % Joint 1 RV
    CobraConfig.ARM_DATA.(cnfgID).SLOW_CALIBRATION_TABLE.Joint1_rev_regions.Text =  array2text(motorMapJ1Rv.(fldID)(1,:));
    CobraConfig.ARM_DATA.(cnfgID).SLOW_CALIBRATION_TABLE.Joint1_rev_stepsizes.Text = array2text(motorMapJ1Rv.(fldID)(2,:));
    % Joint 2 RV
    CobraConfig.ARM_DATA.(cnfgID).SLOW_CALIBRATION_TABLE.Joint2_rev_regions.Text = array2text(motorMapJ2Rv.(fldID)(1,:));
    CobraConfig.ARM_DATA.(cnfgID).SLOW_CALIBRATION_TABLE.Joint2_rev_stepsizes.Text = array2text(motorMapJ2Rv.(fldID)(2,:));

end
 

cobraCfg2xml(CobraConfig, xmlEnd);