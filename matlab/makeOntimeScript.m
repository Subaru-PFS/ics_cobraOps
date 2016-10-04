function output = makeOntimeScript(filename,pids,invertedThetas,invertedPhis,Texp,Tstrobe,ontimes,thetaRefSteps,thetaRefOt,phiRefSteps,phiRefOt,stage,cntrStrkSteps,cntrStrkOt,romSteps,initSteps,initOt)

%% Ontime Tuning Script Generator
% Open file
tfile = fopen(filename,'w');

% Check which direction we are moving
    if romSteps > 0
        drct = 'fw';
    else
        drct = 'rv';
    end

% Create directory and initialize LED
fprintf(tfile,['{ Init camera and LED }\n'...
'cmd_init_Digital_IO\n'...
'cmd_turnLed_Off\n'...
'\n'...
'{Init dir}\n'...
['cmd_createTestDirectoryWithSuffix S' num2str(stage) drct '_onTimeTuning\n']...
'\n'...
'{Set exposure time}\n'...
['cmd_remote_setCameraExposureTime ' num2str(Texp) ', 100\n']...
sprintf('cmd_loadDark_Image  %dsDarkImage_64.fits\n',Texp)...
'cmd_loadBias_Image  RemoteBiasImage_64.fits\n'...
'cmd_enableBiasSub  1\n'...
'cmd_enableDarkSub  1\n'...
'\n'...
'{Init positioners}\n'...
'cmd_zeroPositionerIds\n']);

c = 0;
for pid = pids
    fprintf(tfile,'cmd_setPositionerIds2  %d, 1, %d, 1.4, 2.8, 1.2, 2.4\n',c,pid);
    
    % Check if theta is inverted
    if ~isempty(find(invertedThetas==pid))
        invTheta = 1;
    else
        invTheta = 0;
    end
    
    % Check if phi is inverted
    if ~isempty(find(invertedPhis==pid))
        invPhi = 1;
    else
        invPhi = 0;
    end
    
    fprintf(tfile,'cmd_switchMotorPolarity  1,%d,%d,%d\n',pid,invTheta,invPhi);
    c = c+1;
end

fprintf(tfile,[
'\n'...
'{ Move Theta and Phi to reference positions}\n']);

% Move Phis to reference position
for pid = pids
    fprintf(tfile,['cmd_moveMotor_DurationInterval 1,' num2str(pid) ',2,' num2str(phiRefSteps) ',' num2str(phiRefOt) ',2.5\n']);
end
% fprintf(tfile,['cmd_waitSeconds ' num2str(2.5*abs(phiRefSteps)/1000+1) '\n']);

% Move thetas to initial position
for pid = pids
    fprintf(tfile,['cmd_moveMotor_DurationInterval 1,' num2str(pid) ',1,' num2str(thetaRefSteps) ',' num2str(thetaRefOt) ',2.5\n']);
end
fprintf(tfile,['cmd_waitSeconds ' num2str(2.5*abs(thetaRefSteps)/1000+3) '\n']);


fprintf(tfile,[
'\n'...% Do a reverse streak to find centers with
'{ Take streak exposure for finding stage of interest centers}\n'...
'cmd_getImageStart_NoWait\n'...
'cmd_waitSeconds 0.5\n']);
for pid=pids
    fprintf(tfile,['cmd_moveMotor_DurationInterval 1,' num2str(pid) ',' num2str(stage) ',' num2str(cntrStrkSteps) ',' num2str(cntrStrkOt) ',2.5\n']);
end
fprintf(tfile,[
['cmd_waitSeconds ' num2str(2.5*abs(cntrStrkSteps)/1000+1) '\n']...
'cmd_getImageDone_Wait\n'...
'cmd_saveCurrent_RemoteImage_WithPrefix cntrStrk_\n'...
'cmd_clearOpenCV_Images\n'...
'cmd_waitSeconds  1\n\n']);

%% USING STROBES
% for ot=ontimes
%     fprintf(tfile,[strcat('\n{ ',num2str(100*ot),'MS ONTIME STREAKS}\n')...
%         'cmd_getImageStart_NoWait\n'...
%         'cmd_waitSeconds  2\n'...
%         'cmd_turnLed_On\n'...
%         'cmd_turnLed_Off\n']);
%     
% %             ['cmd_waitSeconds ' num2str(Tstrobe) '\n']...
% 
%     % Run ROM test for each positionger
%     for pid=1:9
%         fprintf(tfile,['cmd_moveMotor_DurationInterval 1,' num2str(pid) ',' num2str(stage) ',' num2str(romSteps) ',' num2str(ot) ',2.5\n']);
%     end
%     
%     fprintf(tfile,[
%     ['cmd_waitSeconds ' num2str(abs(romSteps)*2.5/1000+1) '\n']...
%     'cmd_turnLed_On\n'...
%     'cmd_turnLed_Off\n'...
%     'cmd_getImageDone_Wait\n']);
% 
%     if romSteps > 0
%         drct = 'fw';
%     else
%         drct = 'rv';
%     end
%     fprintf(tfile,['cmd_saveCurrent_RemoteImage_WithPrefix S' num2str(stage) '_ontime_' num2str(ot*1000) 'us_' drct '_\n']);
% 
%     % Return each positioner to initial position
%     for pid=1:9
%         fprintf(tfile,['cmd_moveMotor_DurationInterval 1,' num2str(pid) ',' num2str(stage) ',' num2str(initSteps) ',' num2str(initOt) ',2.5\n']);
%     end
%     fprintf(tfile, ['cmd_waitSeconds ' num2str(abs(initSteps)*2.5/1000+1) '\n']);
%     
%     
% end


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


%% USING STREAKS WITH PAUSES
Nmoves = 20;%PULL ME OUT EVENTUALLY!
smlMoveSteps = round(romSteps/Nmoves); 

for ot=ontimes
    fprintf(tfile,[strcat('\n{ ',num2str(100*ot),'MS ONTIME STREAKS}\n')...
        'cmd_getImageStart_NoWait\n'...
        'cmd_waitSeconds  2\n'...
        ['cmd_Loop_Cp_StartCmd  1, ' num2str(Nmoves) '\n']]);

    % Run ROM test for each positionger
    for pid=pids
        fprintf(tfile,['cmd_moveMotor_DurationInterval 1,' num2str(pid) ',' num2str(stage) ',' num2str(smlMoveSteps) ',' num2str(ot) ',2.5\n']);
    end
    
    fprintf(tfile,[
    ['cmd_waitSeconds ' num2str(abs(smlMoveSteps)*2.5/1000) '\n']...
    'cmd_waitSeconds .1\n'...
    'cmd_Loop_Cp_EndCmd  1\n'...
    'cmd_getImageDone_Wait\n'...
    'cmd_waitSeconds  2\n']);

    fprintf(tfile,'cmd_saveCurrent_RemoteImage_WithPrefix S%d_ontime_%03.0f_us_%s_\n',stage,ot*1000,drct);

    % Return each positioner to initial position
    for pid=pids
        fprintf(tfile,['cmd_moveMotor_DurationInterval 1,' num2str(pid) ',' num2str(stage) ',' num2str(initSteps) ',' num2str(initOt) ',2.5\n']);
    end
    fprintf(tfile, [
        ['cmd_waitSeconds ' num2str(abs(initSteps)*2.5/1000+1) '\n']...
        'cmd_waitSeconds  2\n']);
    
    
end


fclose(tfile);
end

