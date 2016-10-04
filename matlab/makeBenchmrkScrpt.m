%% Motor Mapping MSIM Script Creation

%% INPUTS

% Choose the Config XML to begin with
StartCfg = loadCfgXml;

% Which positioners will be mapped?
mapPosIds = 1:9;

%% EXECUTION

% Create new motor mapping script
[FileName,PathName] = uiputfile('benchmark.lst','Save new motor mapping script as');
tfile = fopen(fullfile(PathName,FileName),'w');

% Get number of ARM_DATA stuctures
Nads = length(find(1-cellfun(@isempty,strfind(fields(StartCfg.ARM_DATA),'ARM_DATA'))));

refpos = str2num(StartCfg.ARM_DATA.refpos.Text);

fprintf(tfile,['{ Init camera and LED }\n'...
'{ Daily Benchmark Test}\n'...
'cmd_init_Digital_IO\n'...
'cmd_turnLed_Off\n'...
'\n'...
'cmd_createTestDirectoryWithSuffix msimMaps\n'...
'\n'...
'{ ==Camera Setup== \n'...
' NOTE: You must run or load dark and bias images prior to this script }\n'...
'cmd_remoteEnableCameraRecovery 1,100\n'...
'cmd_remote_setCameraExposureTime  0.002, 100\n'...
'cmd_selectCentroiding_Algorithm  6\n'...
'cmd_enableWindowed_Centroid  0\n'...
'cmd_setThreshold_NSigma  5\n'...
'cmd_loadDark_Image  2msDarkImage_64.fits\n'...
'cmd_loadBias_Image  RemoteBiasImage_64.fits\n'...
'cmd_enableBiasSub  1\n'...
'cmd_enableDarkSub  1\n'...
'cmd_clearFiducialTable\n'])

for ii=1:length(refpos)
    fprintf(tfile,'cmd_setHornMethodFiducialCoordinate %4.4f, %4.4f, 10\n',real(refpos(ii)),imag(refpos(ii)));
end

% STREAKS
for ii=1:nads
    
keyboard







% fprintf(tfile,['{ MOTOR MAPPING SCRIPT'...
% '======================='...
% ''...
% ' Init camera and LED }'...
% 'cmd_init_Digital_IO'...
% 'cmd_turnLed_Off'...
% ''...
% 'cmd_createTestDirectoryWithSuffix msimMaps'...
% ''...
% '{ ==Camera Setup== '...
% ' NOTE: You must run or load dark and bias images prior to this script }'...
% 'cmd_remoteEnableCameraRecovery 1,100'...
% 'cmd_remote_setCameraExposureTime  0.002, 100'...
% 'cmd_selectCentroiding_Algorithm  6'...
% 'cmd_enableWindowed_Centroid  0'...
% 'cmd_setThreshold_NSigma  5'...
% 'cmd_loadDark_Image  2msDarkImage_64.fits'...
% 'cmd_loadBias_Image  RemoteBiasImage_64.fits'...
% 'cmd_enableBiasSub  1'...
% 'cmd_enableDarkSub  1'...
% 'cmd_clearFiducialTable'...
% 'cmd_setHornMethodFiducialCoordinate %4.4f, %4.4f, 10'...
% 'cmd_setHornMethodFiducialCoordinate %4.4f, %4.4f, 10'...
% 'cmd_setHornMethodFiducialCoordinate %4.4f, %4.4f, 10'],fidx,fidy);

% 'cmd_setHornMethodFiducialCoordinate 1644.1081, 1905.1558, 10'...
% 'cmd_setHornMethodFiducialCoordinate 94.7736, 405.1883, 10'...
% 'cmd_setHornMethodFiducialCoordinate 1257.2365, 1530.2563, 10'...
% 'cmd_enableHornMethod'...
% ''...
% '{ Turn LED on}'...
% 'cmd_turnLed_On'...
% ''...
% '{Init positioners}'...
% 'cmd_zeroPositionerIds'...
% 'cmd_setPositionerIds2  0, 1, 1, 1.4, 2.8, 1.2, 2.4'...
% 'cmd_setPositionerIds2  1, 1, 2, 1.4, 2.8, 1.2, 2.4'...
% 'cmd_setPositionerIds2  2, 1, 3, 1.4, 2.8, 1.2, 2.4'...
% 'cmd_setPositionerIds2  3, 1, 4, 1.4, 2.8, 1.2, 2.4'...
% 'cmd_setPositionerIds2  4, 1, 5, 1.4, 2.8, 1.2, 2.4'...
% 'cmd_setPositionerIds2  5, 1, 6, 1.4, 2.8, 1.2, 2.4'...
% 'cmd_setPositionerIds2  6, 1, 7, 1.4, 2.8, 1.2, 2.4'...
% 'cmd_setPositionerIds2  7, 1, 8, 1.4, 2.8, 1.2, 2.4'...
% 'cmd_setPositionerIds2  8, 1, 9, 1.4, 2.8, 1.2, 2.4'...
% 'cmd_switchMotorPolarity  1,1,0,0'...
% 'cmd_switchMotorPolarity  1,2,0,0'...
% 'cmd_switchMotorPolarity  1,3,0,0'...
% 'cmd_switchMotorPolarity  1,4,0,0'...
% 'cmd_switchMotorPolarity  1,5,0,0'...
% 'cmd_switchMotorPolarity  1,6,0,0'...
% 'cmd_switchMotorPolarity  1,7,0,0'...
% 'cmd_switchMotorPolarity  1,8,0,0'...
% 'cmd_switchMotorPolarity  1,9,0,0'...
% ''...
% '{ Move Theta and Phi to home positions}'...
% 'cmd_moveMotor_DurationInterval 1,1,2,-4000,0.15,2.5'...
% 'cmd_moveMotor_DurationInterval 1,2,2,-4000,0.15,2.5'...
% 'cmd_moveMotor_DurationInterval 1,3,2,-4000,0.15,2.5'...
% 'cmd_moveMotor_DurationInterval 1,4,2,-4000,0.15,2.5'...
% 'cmd_moveMotor_DurationInterval 1,5,2,-4000,0.15,2.5'...
% 'cmd_moveMotor_DurationInterval 1,6,2,-4000,0.15,2.5'...
% 'cmd_moveMotor_DurationInterval 1,7,2,-4000,0.15,2.5'...
% 'cmd_moveMotor_DurationInterval 1,8,2,-4000,0.15,2.5'...
% 'cmd_moveMotor_DurationInterval 1,9,2,-4000,0.15,2.5'...
% 'cmd_moveMotor_DurationInterval 1,1,1,-7000,0.3,2.5'...
% 'cmd_moveMotor_DurationInterval 1,2,1,-7000,0.3,2.5'...
% 'cmd_moveMotor_DurationInterval 1,3,1,-7000,0.3,2.5'...
% 'cmd_moveMotor_DurationInterval 1,4,1,-7000,0.3,2.5'...
% 'cmd_moveMotor_DurationInterval 1,5,1,-7000,0.3,2.5'...
% 'cmd_moveMotor_DurationInterval 1,6,1,-7000,0.3,2.5'...
% 'cmd_moveMotor_DurationInterval 1,7,1,-7000,0.3,2.5'...
% 'cmd_moveMotor_DurationInterval 1,8,1,-7000,0.3,2.5'...
% 'cmd_moveMotor_DurationInterval 1,9,1,-7000,0.3,2.5'...
% 'cmd_waitSeconds 20'...
% ''...
% '{ Move phi out slightly for better joint calculations }'...
% 'cmd_moveMotor_DurationInterval 1,1,2,250,0.06,2.5'...
% 'cmd_moveMotor_DurationInterval 1,2,2,250,0.06,2.5'...
% 'cmd_moveMotor_DurationInterval 1,3,2,250,0.06,2.5'...
% 'cmd_moveMotor_DurationInterval 1,4,2,250,0.06,2.5'...
% 'cmd_moveMotor_DurationInterval 1,5,2,250,0.06,2.5'...
% 'cmd_moveMotor_DurationInterval 1,6,2,250,0.06,2.5'...
% 'cmd_moveMotor_DurationInterval 1,7,2,250,0.07,2.5'...
% 'cmd_moveMotor_DurationInterval 1,8,2,250,0.06,2.5'...
% 'cmd_moveMotor_DurationInterval 1,9,2,250,0.06,2.5'...
% 'cmd_waitSeconds 1'...
% ''...
% '{ ==Setup directory== }'...
% 'cmd_setCentroidsArmLogPrefix    thetaFwMap'...
% '{ Take reference image}'...
% 'cmd_getImageStart_NoWait'...
% 'cmd_getImageDone_Wait'...
% 'cmd_saveCurrent_RemoteImage_WithPrefix thetaFwMap'...
% 'cmd_findCentroids'...
% 'cmd_LogCentroids'...
% 'cmd_clearOpenCV_Images'...
% 'cmd_setCurrentPos 1,1'...
% 'cmd_setCurrentPos 1,2'...
% 'cmd_setCurrentPos 1,3'...
% 'cmd_setCurrentPos 1,4'...
% 'cmd_setCurrentPos 1,5'...
% 'cmd_setCurrentPos 1,6'...
% 'cmd_setCurrentPos 1,7'...
% 'cmd_setCurrentPos 1,8'...
% 'cmd_setCurrentPos 1,9'...
% 'cmd_getCentroidsArmInverseKin  '...
% 'cmd_waitSeconds  1'...
% ''...
% '{ Move theta full ROM FWD in discrete steps }'...
% 'cmd_Loop_Cp_StartCmd  1, 200'...
% ''...
% 'cmd_moveMotor_DurationInterval 1,1,1,60,0.12,2.5'...
% 'cmd_moveMotor_DurationInterval 1,2,1,60,0.06,2.5'...
% 'cmd_moveMotor_DurationInterval 1,3,1,60,0.12,2.5'...
% 'cmd_moveMotor_DurationInterval 1,4,1,60,0.06,2.5'...
% 'cmd_moveMotor_DurationInterval 1,5,1,60,0.08,2.5'...
% 'cmd_moveMotor_DurationInterval 1,6,1,60,0.06,2.5'...
% 'cmd_moveMotor_DurationInterval 1,7,1,60,0.14,2.5'...
% 'cmd_moveMotor_DurationInterval 1,8,1,60,0.08,2.5'...
% 'cmd_moveMotor_DurationInterval 1,9,1,60,0.1,2.5'...
% 'cmd_waitSeconds 0.5'...
% ''...
% '{ Take reference image}'...
% 'cmd_getImageStart_NoWait'...
% 'cmd_getImageDone_Wait'...
% 'cmd_saveCurrent_RemoteImage_WithPrefix thetaFwMap'...
% 'cmd_findCentroids'...
% 'cmd_LogCentroids'...
% 'cmd_clearOpenCV_Images'...
% 'cmd_setCurrentPos 1,1'...
% 'cmd_setCurrentPos 1,2'...
% 'cmd_setCurrentPos 1,3'...
% 'cmd_setCurrentPos 1,4'...
% 'cmd_setCurrentPos 1,5'...
% 'cmd_setCurrentPos 1,6'...
% 'cmd_setCurrentPos 1,7'...
% 'cmd_setCurrentPos 1,8'...
% 'cmd_setCurrentPos 1,9'...
% 'cmd_getCentroidsArmInverseKin  '...
% 'cmd_Loop_Cp_EndCmd  1'...
% ''...
% '{ Move Theta to anti-home positions}'...
% 'cmd_moveMotor_DurationInterval 1,1,1,7000,0.3,2.5'...
% 'cmd_moveMotor_DurationInterval 1,2,1,7000,0.3,2.5'...
% 'cmd_moveMotor_DurationInterval 1,3,1,7000,0.3,2.5'...
% 'cmd_moveMotor_DurationInterval 1,4,1,7000,0.3,2.5'...
% 'cmd_moveMotor_DurationInterval 1,5,1,7000,0.3,2.5'...
% 'cmd_moveMotor_DurationInterval 1,6,1,7000,0.3,2.5'...
% 'cmd_moveMotor_DurationInterval 1,7,1,7000,0.3,2.5'...
% 'cmd_moveMotor_DurationInterval 1,8,1,7000,0.3,2.5'...
% 'cmd_moveMotor_DurationInterval 1,9,1,7000,0.3,2.5'...
% 'cmd_waitSeconds 20'...
% ''...
% '{ ==Setup directory== }'...
% 'cmd_setCentroidsArmLogPrefix    thetaRvMap'...
% '{ Take reference image}'...
% 'cmd_getImageStart_NoWait'...
% 'cmd_getImageDone_Wait'...
% 'cmd_saveCurrent_RemoteImage_WithPrefix thetaRvMap'...
% 'cmd_findCentroids'...
% 'cmd_LogCentroids'...
% 'cmd_clearOpenCV_Images'...
% 'cmd_setCurrentPos 1,1'...
% 'cmd_setCurrentPos 1,2'...
% 'cmd_setCurrentPos 1,3'...
% 'cmd_setCurrentPos 1,4'...
% 'cmd_setCurrentPos 1,5'...
% 'cmd_setCurrentPos 1,6'...
% 'cmd_setCurrentPos 1,7'...
% 'cmd_setCurrentPos 1,8'...
% 'cmd_setCurrentPos 1,9'...
% 'cmd_getCentroidsArmInverseKin  '...
% 'cmd_waitSeconds  1'...
% ''...
% '{ Move theta full ROM REV in discrete steps }'...
% 'cmd_Loop_Cp_StartCmd  2, 200'...
% ''...
% 'cmd_moveMotor_DurationInterval 1,1,1,-60,0.1,2.5'...
% 'cmd_moveMotor_DurationInterval 1,2,1,-60,0.08,2.5'...
% 'cmd_moveMotor_DurationInterval 1,3,1,-60,0.1,2.5'...
% 'cmd_moveMotor_DurationInterval 1,4,1,-60,0.08,2.5'...
% 'cmd_moveMotor_DurationInterval 1,5,1,-60,0.08,2.5'...
% 'cmd_moveMotor_DurationInterval 1,6,1,-60,0.08,2.5'...
% 'cmd_moveMotor_DurationInterval 1,7,1,-60,0.14,2.5'...
% 'cmd_moveMotor_DurationInterval 1,8,1,-60,0.08,2.5'...
% 'cmd_moveMotor_DurationInterval 1,9,1,-60,0.08,2.5'...
% 'cmd_waitSeconds 0.5'...
% ''...
% '{ Take reference image}'...
% 'cmd_getImageStart_NoWait'...
% 'cmd_getImageDone_Wait'...
% 'cmd_saveCurrent_RemoteImage_WithPrefix thetaRvMap'...
% 'cmd_findCentroids'...
% 'cmd_LogCentroids'...
% 'cmd_clearOpenCV_Images'...
% 'cmd_setCurrentPos 1,1'...
% 'cmd_setCurrentPos 1,2'...
% 'cmd_setCurrentPos 1,3'...
% 'cmd_setCurrentPos 1,4'...
% 'cmd_setCurrentPos 1,5'...
% 'cmd_setCurrentPos 1,6'...
% 'cmd_setCurrentPos 1,7'...
% 'cmd_setCurrentPos 1,8'...
% 'cmd_setCurrentPos 1,9'...
% 'cmd_getCentroidsArmInverseKin  '...
% 'cmd_Loop_Cp_EndCmd  2'...
% ''...
% '{ Move Theta and Phi to home positions}'...
% 'cmd_moveMotor_DurationInterval 1,1,2,-4000,0.15,2.5'...
% 'cmd_moveMotor_DurationInterval 1,2,2,-4000,0.15,2.5'...
% 'cmd_moveMotor_DurationInterval 1,3,2,-4000,0.15,2.5'...
% 'cmd_moveMotor_DurationInterval 1,4,2,-4000,0.15,2.5'...
% 'cmd_moveMotor_DurationInterval 1,5,2,-4000,0.15,2.5'...
% 'cmd_moveMotor_DurationInterval 1,6,2,-4000,0.15,2.5'...
% 'cmd_moveMotor_DurationInterval 1,7,2,-4000,0.15,2.5'...
% 'cmd_moveMotor_DurationInterval 1,8,2,-4000,0.15,2.5'...
% 'cmd_moveMotor_DurationInterval 1,9,2,-4000,0.15,2.5'...
% 'cmd_moveMotor_DurationInterval 1,1,1,-7000,0.3,2.5'...
% 'cmd_moveMotor_DurationInterval 1,2,1,-7000,0.3,2.5'...
% 'cmd_moveMotor_DurationInterval 1,3,1,-7000,0.3,2.5'...
% 'cmd_moveMotor_DurationInterval 1,4,1,-7000,0.3,2.5'...
% 'cmd_moveMotor_DurationInterval 1,5,1,-7000,0.3,2.5'...
% 'cmd_moveMotor_DurationInterval 1,6,1,-7000,0.3,2.5'...
% 'cmd_moveMotor_DurationInterval 1,7,1,-7000,0.3,2.5'...
% 'cmd_moveMotor_DurationInterval 1,8,1,-7000,0.3,2.5'...
% 'cmd_moveMotor_DurationInterval 1,9,1,-7000,0.3,2.5'...
% 'cmd_waitSeconds 20'...
% ''...
% '{ ==Setup directory== }'...
% 'cmd_setCentroidsArmLogPrefix    PhiFwMap'...
% '{ Take reference image}'...
% 'cmd_getImageStart_NoWait'...
% 'cmd_getImageDone_Wait'...
% 'cmd_saveCurrent_RemoteImage_WithPrefix PhiFwMap'...
% 'cmd_findCentroids'...
% 'cmd_LogCentroids'...
% 'cmd_clearOpenCV_Images'...
% 'cmd_setCurrentPos 1,1'...
% 'cmd_setCurrentPos 1,2'...
% 'cmd_setCurrentPos 1,3'...
% 'cmd_setCurrentPos 1,4'...
% 'cmd_setCurrentPos 1,5'...
% 'cmd_setCurrentPos 1,6'...
% 'cmd_setCurrentPos 1,7'...
% 'cmd_setCurrentPos 1,8'...
% 'cmd_setCurrentPos 1,9'...
% 'cmd_getCentroidsArmInverseKin  '...
% 'cmd_waitSeconds  1'...
% ''...
% '{ Move theta and phi full ROM FWD in discrete steps }'...
% 'cmd_Loop_Cp_StartCmd  3, 150'...
% ''...
% 'cmd_moveMotor_DurationInterval 1,1,2,30,0.06,2.5'...
% 'cmd_moveMotor_DurationInterval 1,2,2,30,0.06,2.5'...
% 'cmd_moveMotor_DurationInterval 1,3,2,30,0.06,2.5'...
% 'cmd_moveMotor_DurationInterval 1,4,2,30,0.06,2.5'...
% 'cmd_moveMotor_DurationInterval 1,5,2,30,0.06,2.5'...
% 'cmd_moveMotor_DurationInterval 1,6,2,30,0.06,2.5'...
% 'cmd_moveMotor_DurationInterval 1,7,2,30,0.07,2.5'...
% 'cmd_moveMotor_DurationInterval 1,8,2,30,0.06,2.5'...
% 'cmd_moveMotor_DurationInterval 1,9,2,30,0.06,2.5'...
% 'cmd_waitSeconds 1'...
% ''...
% '{ Take reference image}'...
% 'cmd_getImageStart_NoWait'...
% 'cmd_getImageDone_Wait'...
% 'cmd_saveCurrent_RemoteImage_WithPrefix PhiFwMap'...
% 'cmd_findCentroids'...
% 'cmd_LogCentroids'...
% 'cmd_clearOpenCV_Images'...
% 'cmd_setCurrentPos 1,1'...
% 'cmd_setCurrentPos 1,2'...
% 'cmd_setCurrentPos 1,3'...
% 'cmd_setCurrentPos 1,4'...
% 'cmd_setCurrentPos 1,5'...
% 'cmd_setCurrentPos 1,6'...
% 'cmd_setCurrentPos 1,7'...
% 'cmd_setCurrentPos 1,8'...
% 'cmd_setCurrentPos 1,9'...
% 'cmd_getCentroidsArmInverseKin  '...
% 'cmd_waitSeconds  1'...
% 'cmd_Loop_Cp_EndCmd  3'...
% ''...
% '{ Move Phi to antihome position}'...
% 'cmd_moveMotor_DurationInterval 1,1,2,4000,0.15,2.5'...
% 'cmd_moveMotor_DurationInterval 1,2,2,4000,0.15,2.5'...
% 'cmd_moveMotor_DurationInterval 1,3,2,4000,0.15,2.5'...
% 'cmd_moveMotor_DurationInterval 1,4,2,4000,0.15,2.5'...
% 'cmd_moveMotor_DurationInterval 1,5,2,4000,0.15,2.5'...
% 'cmd_moveMotor_DurationInterval 1,6,2,4000,0.15,2.5'...
% 'cmd_moveMotor_DurationInterval 1,7,2,4000,0.15,2.5'...
% 'cmd_moveMotor_DurationInterval 1,8,2,4000,0.15,2.5'...
% 'cmd_moveMotor_DurationInterval 1,9,2,4000,0.15,2.5'...
% 'cmd_waitSeconds 13'...
% ''...
% '{ ==Setup directory== }'...
% 'cmd_setCentroidsArmLogPrefix    PhiRvMap'...
% '{ Take reference image}'...
% 'cmd_getImageStart_NoWait'...
% 'cmd_getImageDone_Wait'...
% 'cmd_saveCurrent_RemoteImage_WithPrefix PhiRvMap'...
% 'cmd_findCentroids'...
% 'cmd_LogCentroids'...
% 'cmd_clearOpenCV_Images'...
% 'cmd_setCurrentPos 1,1'...
% 'cmd_setCurrentPos 1,2'...
% 'cmd_setCurrentPos 1,3'...
% 'cmd_setCurrentPos 1,4'...
% 'cmd_setCurrentPos 1,5'...
% 'cmd_setCurrentPos 1,6'...
% 'cmd_setCurrentPos 1,7'...
% 'cmd_setCurrentPos 1,8'...
% 'cmd_setCurrentPos 1,9'...
% 'cmd_getCentroidsArmInverseKin  '...
% 'cmd_waitSeconds  1'...
% ''...
% '{ Move theta and phi full ROM REV in discrete steps }'...
% 'cmd_Loop_Cp_StartCmd  4, 150'...
% ''...
% 'cmd_moveMotor_DurationInterval 1,1,2,-30,0.05,2.5'...
% 'cmd_moveMotor_DurationInterval 1,2,2,-30,0.05,2.5'...
% 'cmd_moveMotor_DurationInterval 1,3,2,-30,0.05,2.5'...
% 'cmd_moveMotor_DurationInterval 1,4,2,-30,0.06,2.5'...
% 'cmd_moveMotor_DurationInterval 1,5,2,-30,0.07,2.5'...
% 'cmd_moveMotor_DurationInterval 1,6,2,-30,0.05,2.5'...
% 'cmd_moveMotor_DurationInterval 1,7,2,-30,0.06,2.5'...
% 'cmd_moveMotor_DurationInterval 1,8,2,-30,0.05,2.5'...
% 'cmd_moveMotor_DurationInterval 1,9,2,-30,0.05,2.5'...
% 'cmd_waitSeconds 1'...
% ''...
% '{ Take reference image}'...
% 'cmd_getImageStart_NoWait'...
% 'cmd_getImageDone_Wait'...
% 'cmd_saveCurrent_RemoteImage_WithPrefix PhiRvMap'...
% 'cmd_findCentroids'...
% 'cmd_LogCentroids'...
% 'cmd_clearOpenCV_Images'...
% 'cmd_setCurrentPos 1,1'...
% 'cmd_setCurrentPos 1,2'...
% 'cmd_setCurrentPos 1,3'...
% 'cmd_setCurrentPos 1,4'...
% 'cmd_setCurrentPos 1,5'...
% 'cmd_setCurrentPos 1,6'...
% 'cmd_setCurrentPos 1,7'...
% 'cmd_setCurrentPos 1,8'...
% 'cmd_setCurrentPos 1,9'...
% 'cmd_getCentroidsArmInverseKin  '...
% 'cmd_waitSeconds  1'...
% 'cmd_Loop_Cp_EndCmd  4'...
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% % Create directory and initialize LED
% fprintf(tfile,['{ Init camera and LED }\n'...
% 'cmd_init_Digital_IO\n'...
% 'cmd_turnLed_Off\n'...
% '\n'...
% '{Init dir}\n'...
% 'cmd_createTestDirectoryWithSuffix onTimeTuning\n'...
% '\n'...
% '{Set exposure time}\n'...
% ['cmd_remote_setCameraExposureTime ' num2str(Texp) ', 100\n']...
% sprintf('cmd_loadDark_Image  %dsDarkImage_64.fits\n',Texp)...
% 'cmd_loadBias_Image  RemoteBiasImage_64.fits\n'...
% 'cmd_enableBiasSub  1\n'...
% 'cmd_enableDarkSub  1\n'...
% '\n'...
% '{Init positioners}\n'...
% 'cmd_zeroPositionerIds\n']);
% 
% c = 0;
% for pid = pids
%     fprintf(tfile,'cmd_setPositionerIds2  %d, 1, %d, 1.4, 2.8, 1.2, 2.4\n',c,pid);
%     
%     % Check if theta is inverted
%     if ~isempty(find(invertedThetas==pid))
%         invTheta = 1;
%     else
%         invTheta = 0;
%     end
%     
%     % Check if phi is inverted
%     if ~isempty(find(invertedPhis==pid))
%         invPhi = 1;
%     else
%         invPhi = 0;
%     end
%     
%     fprintf(tfile,'cmd_switchMotorPolarity  1,%d,%d,%d\n',pid,invTheta,invPhi);
%     c = c+1;
% end
% 
% fprintf(tfile,[
% '\n'...
% '{ Move Theta and Phi to reference positions}\n']);
% 
% % Move Phis to reference position
% for pid = pids
%     fprintf(tfile,['cmd_moveMotor_DurationInterval 1,' num2str(pid) ',2,' num2str(phiRefSteps) ',' num2str(phiRefOt) ',2.5\n']);
% end
% % fprintf(tfile,['cmd_waitSeconds ' num2str(2.5*abs(phiRefSteps)/1000+1) '\n']);
% 
% % Move thetas to initial position
% for pid = pids
%     fprintf(tfile,['cmd_moveMotor_DurationInterval 1,' num2str(pid) ',1,' num2str(thetaRefSteps) ',' num2str(thetaRefOt) ',2.5\n']);
% end
% fprintf(tfile,['cmd_waitSeconds ' num2str(2.5*abs(thetaRefSteps)/1000+3) '\n']);
% 
% 
% fprintf(tfile,[
% '\n'...% Do a reverse streak to find centers with
% '{ Take streak exposure for finding stage of interest centers}\n'...
% 'cmd_getImageStart_NoWait\n'...
% 'cmd_waitSeconds 0.5\n']);
% for pid=pids
%     fprintf(tfile,['cmd_moveMotor_DurationInterval 1,' num2str(pid) ',' num2str(stage) ',' num2str(cntrStrkSteps) ',' num2str(cntrStrkOt) ',2.5\n']);
% end
% fprintf(tfile,[
% ['cmd_waitSeconds ' num2str(2.5*abs(cntrStrkSteps)/1000+1) '\n']...
% 'cmd_getImageDone_Wait\n'...
% 'cmd_saveCurrent_RemoteImage_WithPrefix cntrStrk_\n'...
% 'cmd_clearOpenCV_Images\n'...
% 'cmd_waitSeconds  1\n\n']);
% 
% %% USING STROBES
% % for ot=ontimes
% %     fprintf(tfile,[strcat('\n{ ',num2str(100*ot),'MS ONTIME STREAKS}\n')...
% %         'cmd_getImageStart_NoWait\n'...
% %         'cmd_waitSeconds  2\n'...
% %         'cmd_turnLed_On\n'...
% %         'cmd_turnLed_Off\n']);
% %     
% % %             ['cmd_waitSeconds ' num2str(Tstrobe) '\n']...
% % 
% %     % Run ROM test for each positionger
% %     for pid=1:9
% %         fprintf(tfile,['cmd_moveMotor_DurationInterval 1,' num2str(pid) ',' num2str(stage) ',' num2str(romSteps) ',' num2str(ot) ',2.5\n']);
% %     end
% %     
% %     fprintf(tfile,[
% %     ['cmd_waitSeconds ' num2str(abs(romSteps)*2.5/1000+1) '\n']...
% %     'cmd_turnLed_On\n'...
% %     'cmd_turnLed_Off\n'...
% %     'cmd_getImageDone_Wait\n']);
% % 
% %     if romSteps > 0
% %         drct = 'fw';
% %     else
% %         drct = 'rv';
% %     end
% %     fprintf(tfile,['cmd_saveCurrent_RemoteImage_WithPrefix S' num2str(stage) '_ontime_' num2str(ot*1000) 'us_' drct '_\n']);
% % 
% %     % Return each positioner to initial position
% %     for pid=1:9
% %         fprintf(tfile,['cmd_moveMotor_DurationInterval 1,' num2str(pid) ',' num2str(stage) ',' num2str(initSteps) ',' num2str(initOt) ',2.5\n']);
% %     end
% %     fprintf(tfile, ['cmd_waitSeconds ' num2str(abs(initSteps)*2.5/1000+1) '\n']);
% %     
% %     
% % end
% 
% 
% %% USING STREAKS
% for ot=ontimes
%     fprintf(tfile,[strcat('\n{ ',num2str(100*ot),'MS ONTIME STREAKS}\n')...
%         'cmd_getImageStart_NoWait\n'...
%         'cmd_turnLed_On\n'...
%         'cmd_waitSeconds  2\n']);
% 
%     % Run ROM test for each positionger
%     for pid=pids
%         fprintf(tfile,['cmd_moveMotor_DurationInterval 1,' num2str(pid) ',' num2str(stage) ',' num2str(romSteps) ',' num2str(ot) ',2.5\n']);
%     end
%     
%     fprintf(tfile,[
%     ['cmd_waitSeconds ' num2str(abs(romSteps)*2.5/1000+1) '\n']...
%     'cmd_turnLed_Off\n'...
%     'cmd_getImageDone_Wait\n'...
%     'cmd_waitSeconds  2\n']);
% 
%     if romSteps > 0
%         drct = 'fw';
%     else
%         drct = 'rv';
%     end
% 
%     fprintf(tfile,'cmd_saveCurrent_RemoteImage_WithPrefix S%d_ontime_%03.0f_us_%s_\n',stage,ot*1000,drct);
% 
%     % Return each positioner to initial position
%     for pid=pids
%         fprintf(tfile,['cmd_moveMotor_DurationInterval 1,' num2str(pid) ',' num2str(stage) ',' num2str(initSteps) ',' num2str(initOt) ',2.5\n']);
%     end
%     fprintf(tfile, [
%         ['cmd_waitSeconds ' num2str(abs(initSteps)*2.5/1000+1) '\n']...
%         'cmd_waitSeconds  2\n']);
%     
%     
% end



fclose(tfile);
% end

