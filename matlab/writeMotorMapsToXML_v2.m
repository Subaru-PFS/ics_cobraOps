%Specifiy XML file
clear all
close all

% Choose xml file to start with
CobraConfig = loadCfgXml;

xmlEnd = fullfile(CobraConfig.cfgPath,['new_v2_' CobraConfig.cfgFile]);

%% load 4 motor maps
mmap_path = fullfile(getenv('HOME'),'Dropbox','PFS_EM','TEST_RESULTS','StreakResults');
load(fullfile(mmap_path,'workspacestage1.mat'),'motorMapFw','motorMapRv');
J1.fw.motorMap = motorMapFw;
J1.rv.motorMap = motorMapRv;
load(fullfile(mmap_path,'workspacestage2.mat'),'motorMapFw','motorMapRv');
J2.fw.motorMap = motorMapFw;
J2.rv.motorMap = motorMapRv;
clear motorMapRv motorMapFw

return

%% FIX AutofindOntimes.m first

J1.fw.OnTimes  = [.07,.07,.13,.07,.18,.25,.07];
J1.rv.OnTimes  = [.13,.09,.13,.07,.18,.30,.09];
J2.fw.OnTimes  = [.05,.05,.05,.05,.05,.05,.05];
J2.rv.OnTimes  = [.05,.05,.05,.05,.05,.05,.05];

J1.fw.IntTimes = [2.5,2.5,2.5,2.5,2.5,2.5,2.5];
J1.rv.IntTimes = [2.5,2.5,2.5,2.5,2.5,2.5,2.5];
J2.fw.IntTimes = [2.5,2.5,2.5,2.5,2.5,2.5,2.5];
J2.rv.IntTimes = [2.5,2.5,2.5,2.5,2.5,2.5,2.5];

for pid=1:7
    
    if pid > 4
        cpid = pid+2;
    else
        cpid = pid;
    end
    
    fldID = sprintf('pId%d',pid);
    cnfgID = sprintf('ARM_DATA_%d',cpid);
    
    CobraConfig = setARMval(CobraConfig,cpid,'Link1_fwd_Duration'      ,J1.fw.OnTimes(pid));
    CobraConfig = setARMval(CobraConfig,cpid,'Link1_fwd_Duration_Slow' ,J1.fw.OnTimes(pid)); 
    CobraConfig = setARMval(CobraConfig,cpid,'Link1_fwd_Intervals'     ,J1.fw.IntTimes(pid));
    CobraConfig = setARMval(CobraConfig,cpid,'Link1_fwd_Intervals_Slow',J1.fw.IntTimes(pid));
    CobraConfig = setARMval(CobraConfig,cpid,'Link1_rev_Duration'      ,J1.rv.OnTimes(pid));
    CobraConfig = setARMval(CobraConfig,cpid,'Link1_rev_Duration_Slow' ,J1.rv.OnTimes(pid));
    CobraConfig = setARMval(CobraConfig,cpid,'Link1_rev_Intervals'     ,J1.rv.IntTimes(pid));
    CobraConfig = setARMval(CobraConfig,cpid,'Link1_rev_Intervals_Slow',J1.rv.IntTimes(pid));

    CobraConfig = setARMval(CobraConfig,cpid,'Link2_fwd_Duration'      ,J2.fw.OnTimes(pid));
    CobraConfig = setARMval(CobraConfig,cpid,'Link2_fwd_Duration_Slow' ,J2.fw.OnTimes(pid));
    CobraConfig = setARMval(CobraConfig,cpid,'Link2_fwd_Intervals'     ,J2.fw.IntTimes(pid));
    CobraConfig = setARMval(CobraConfig,cpid,'Link2_fwd_Intervals_Slow',J2.fw.IntTimes(pid));
    CobraConfig = setARMval(CobraConfig,cpid,'Link2_rev_Duration'      ,J2.rv.OnTimes(pid));
    CobraConfig = setARMval(CobraConfig,cpid,'Link2_rev_Duration_Slow' ,J2.rv.OnTimes(pid));
    CobraConfig = setARMval(CobraConfig,cpid,'Link2_rev_Intervals'     ,J2.rv.IntTimes(pid));
    CobraConfig = setARMval(CobraConfig,cpid,'Link2_rev_Intervals_Slow',J2.rv.IntTimes(pid));

    CobraConfig = setARMval(CobraConfig,cpid,'Joint1_transition_angle' ,400);
    CobraConfig = setARMval(CobraConfig,cpid,'Joint2_transition_angle' ,400);
    
    %%FAST MOTOR MAPS
    %Joint 1 FW
    CobraConfig.ARM_DATA.(cnfgID).FAST_CALIBRATION_TABLE.Joint1_fwd_regions.Text   = array2text(J1.fw.motorMap.(fldID)(1,:));
    CobraConfig.ARM_DATA.(cnfgID).FAST_CALIBRATION_TABLE.Joint1_fwd_stepsizes.Text = array2text(J1.fw.motorMap.(fldID)(2,:));
    % Joint 2 FW
    CobraConfig.ARM_DATA.(cnfgID).FAST_CALIBRATION_TABLE.Joint2_fwd_regions.Text   = array2text(J2.fw.motorMap.(fldID)(1,:));
    CobraConfig.ARM_DATA.(cnfgID).FAST_CALIBRATION_TABLE.Joint2_fwd_stepsizes.Text = array2text(J2.fw.motorMap.(fldID)(2,:));

    % Joint 1 RV
    CobraConfig.ARM_DATA.(cnfgID).FAST_CALIBRATION_TABLE.Joint1_rev_regions.Text   = array2text(J1.rv.motorMap.(fldID)(1,:));
    CobraConfig.ARM_DATA.(cnfgID).FAST_CALIBRATION_TABLE.Joint1_rev_stepsizes.Text = array2text(J1.rv.motorMap.(fldID)(2,:));
    % Joint 2 RV
    CobraConfig.ARM_DATA.(cnfgID).FAST_CALIBRATION_TABLE.Joint2_rev_regions.Text   = array2text(J2.rv.motorMap.(fldID)(1,:));
    CobraConfig.ARM_DATA.(cnfgID).FAST_CALIBRATION_TABLE.Joint2_rev_stepsizes.Text = array2text(J2.rv.motorMap.(fldID)(2,:));


    %% SLOW MOTOR MAPS
    %Joint 1 FW
    CobraConfig.ARM_DATA.(cnfgID).SLOW_CALIBRATION_TABLE.Joint1_fwd_regions.Text   = array2text(J1.fw.motorMap.(fldID)(1,:));
    CobraConfig.ARM_DATA.(cnfgID).SLOW_CALIBRATION_TABLE.Joint1_fwd_stepsizes.Text = array2text(J1.fw.motorMap.(fldID)(2,:));
    % Joint 2 FW
    CobraConfig.ARM_DATA.(cnfgID).SLOW_CALIBRATION_TABLE.Joint2_fwd_regions.Text   = array2text(J2.fw.motorMap.(fldID)(1,:));
    CobraConfig.ARM_DATA.(cnfgID).SLOW_CALIBRATION_TABLE.Joint2_fwd_stepsizes.Text = array2text(J2.fw.motorMap.(fldID)(2,:));

    % Joint 1 RV
    CobraConfig.ARM_DATA.(cnfgID).SLOW_CALIBRATION_TABLE.Joint1_rev_regions.Text   = array2text(J1.rv.motorMap.(fldID)(1,:));
    CobraConfig.ARM_DATA.(cnfgID).SLOW_CALIBRATION_TABLE.Joint1_rev_stepsizes.Text = array2text(J1.rv.motorMap.(fldID)(2,:));
    % Joint 2 RV
    CobraConfig.ARM_DATA.(cnfgID).SLOW_CALIBRATION_TABLE.Joint2_rev_regions.Text   = array2text(J2.rv.motorMap.(fldID)(1,:));
    CobraConfig.ARM_DATA.(cnfgID).SLOW_CALIBRATION_TABLE.Joint2_rev_stepsizes.Text = array2text(J2.rv.motorMap.(fldID)(2,:));

end
 

cobraCfg2xml(CobraConfig, xmlEnd);