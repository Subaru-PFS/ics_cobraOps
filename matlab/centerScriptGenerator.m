function output = ontimeTuning(stage, dir, filename)

Texp = 25;%sec
Tstrobe = 0.5;%sec
phiRefOt = 0.15;
thetaRefOt = .3;

switch stage
    % THETA
    case 1
        ontimes = 0.05:0.02:0.30;
        cntrStrkOt = 0.2;
        initOt = 0.3;
        switch dir
            % FWD
            case 1
                thetaRefSteps = 4000;
                phiRefSteps = 4000;
                cntrStrkSteps = -4000;
                romSteps = 6000;
                initSteps = -4000;
            % REV
            case 2
                thetaRefSteps = -4000;
                phiRefSteps = 4000;
                cntrStrkSteps = 4000;
                romSteps = -6000;
                initSteps = 4000;
        end
    
    % PHI
    case 2
        ontimes = 0.05:0.01:0.15;
        cntrStrkOt = 0.15;
        initOt = 0.15;
        switch dir
            % FWD
            case 1
                thetaRefSteps = 4000;
                phiRefSteps = 4000;
                cntrStrkSteps = -4000;
                romSteps = 3000;
                initSteps = -4000;
            % REV
            case 2
                thetaRefSteps = -4000;
                phiRefSteps = -4000;
                cntrStrkSteps = 4000;
                romSteps = -3000;
                initSteps = 4000;
        end
end

makeOntimeScript(filename,Texp,Tstrobe,ontimes,thetaRefSteps,thetaRefOt,phiRefSteps,phiRefOt,stage,cntrStrkSteps,cntrStrkOt,romSteps,initSteps,initOt);

output = 0;

end
                



function output = makeOntimeScript(filename,Texp,Tstrobe,ontimes,thetaRefSteps,thetaRefOt,phiRefSteps,phiRefOt,stage,cntrStrkSteps,cntrStrkOt,romSteps,initSteps,initOt)

%% Ontime Tuning Script Generator
% Open file
tfile = fopen(filename,'w');

% Create directory and initialize LED
fprintf(tfile,['{ Init camera and LED }\n'...
'cmd_init_Digital_IO\n'...
'cmd_turnLed_Off\n'...
'\n'...
'{Init dir}\n'...
'cmd_createTestDirectory\n'...
'\n'...
'{Set exposure time}\n'...
['cmd_remote_setCameraExposureTime ' num2str(Texp) ', 100\n']...
'\n'...
'\n'...
'{Init positioners}\n'...
'cmd_zeroPositionerIds\n'...
'cmd_setPositionerIds2  0, 1, 1, 1.4, 2.8, 1.2, 2.4\n'...
'cmd_setPositionerIds2  1, 1, 2, 1.4, 2.8, 1.2, 2.4\n'...
'cmd_setPositionerIds2  2, 1, 3, 1.4, 2.8, 1.2, 2.4\n'...
'cmd_setPositionerIds2  3, 1, 4, 1.4, 2.8, 1.2, 2.4\n'...
'cmd_setPositionerIds2  4, 1, 5, 1.4, 2.8, 1.2, 2.4\n'...
'cmd_setPositionerIds2  5, 1, 6, 1.4, 2.8, 1.2, 2.4\n'...
'cmd_setPositionerIds2  6, 1, 7, 1.4, 2.8, 1.2, 2.4\n'...
'cmd_setPositionerIds2  7, 1, 8, 1.4, 2.8, 1.2, 2.4\n'...
'cmd_setPositionerIds2  8, 1, 9, 1.4, 2.8, 1.2, 2.4\n'...
'cmd_switchMotorPolarity  1,1,0,0\n'...
'cmd_switchMotorPolarity  1,2,0,0\n'...
'cmd_switchMotorPolarity  1,3,0,0\n'...
'cmd_switchMotorPolarity  1,4,1,0\n'...
'cmd_switchMotorPolarity  1,5,0,0\n'...
'cmd_switchMotorPolarity  1,6,0,0\n'...
'cmd_switchMotorPolarity  1,7,0,0\n'...
'cmd_switchMotorPolarity  1,8,1,0\n'...
'cmd_switchMotorPolarity  1,9,0,0\n'...
'\n'...
'{ Move Theta and Phi to reference positions}\n']);

% Move Phis to reference position
for pid = [1:4 6:9]
    fprintf(tfile,['cmd_moveMotor_DurationInterval 1,' num2str(pid) ',2,' num2str(phiRefSteps) ',' num2str(phiRefOt) ',2.5\n']);
end
% fprintf(tfile,['cmd_waitSeconds ' num2str(2.5*abs(phiRefSteps)/1000+1) '\n']);

% Move thetas to initial position
for pid = [1:4 6:9]
    fprintf(tfile,['cmd_moveMotor_DurationInterval 1,' num2str(pid) ',1,' num2str(thetaRefSteps) ',' num2str(thetaRefOt) ',2.5\n']);
end
fprintf(tfile,['cmd_waitSeconds ' num2str(2.5*abs(thetaRefSteps)/1000+3) '\n']);

% Take reference image
fprintf(tfile,[
'\n'...
'{ Take reference image}\n'...
'cmd_getImageStart_NoWait\n'...
'cmd_waitSeconds  2\n'...
'cmd_turnLed_On\n'...
'cmd_waitSeconds  ',num2str(Tstrobe),'\n'...
'cmd_turnLed_Off\n'...
'cmd_getImageDone_Wait\n'...
'cmd_saveCurrent_RemoteImage_WithPrefix refPos_\n'...
'cmd_clearOpenCV_Images\n'...
'cmd_turnLed_On\n'...
'cmd_waitSeconds  2\n'...
'\n'...% Do a reverse streak to find centers with
'{ Take streak exposure for finding stage of interest centers}\n'...
'cmd_getImageStart_NoWait\n'...
'cmd_waitSeconds 0.5\n']);
for pid=1:9
    fprintf(tfile,['cmd_moveMotor_DurationInterval 1,' num2str(pid) ',' num2str(stage) ',' num2str(cntrStrkSteps) ',' num2str(cntrStrkOt) ',2.5\n']);
end
fprintf(tfile,[
['cmd_waitSeconds ' num2str(2.5*abs(cntrStrkSteps)/1000+1) '\n']...
'cmd_getImageDone_Wait\n'...
'cmd_saveCurrent_RemoteImage_WithPrefix cntrStrk_\n'...
'cmd_clearOpenCV_Images\n'...
'cmd_waitSeconds  2\n\n']);

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


%% USING STREAKS
for ot=ontimes
    fprintf(tfile,[strcat('\n{ ',num2str(100*ot),'MS ONTIME STREAKS}\n')...
        'cmd_getImageStart_NoWait\n'...
        'cmd_turnLed_On\n'...
        'cmd_waitSeconds  3\n']);

    % Run ROM test for each positionger
    for pid=[1:4 6:9]
        fprintf(tfile,['cmd_moveMotor_DurationInterval 1,' num2str(pid) ',' num2str(stage) ',' num2str(romSteps) ',' num2str(ot) ',2.5\n']);
    end
    
    fprintf(tfile,[
    ['cmd_waitSeconds ' num2str(abs(romSteps)*2.5/1000+1) '\n']...
    'cmd_turnLed_Off\n'...
    'cmd_getImageDone_Wait\n'...
    'cmd_waitSeconds  2\n']);

    if romSteps > 0
        drct = 'fw';
    else
        drct = 'rv';
    end
    fprintf(tfile,['cmd_saveCurrent_RemoteImage_WithPrefix S' num2str(stage) '_ontime_' num2str(ot*1000) 'us_' drct '_\n']);

    % Return each positioner to initial position
    for pid=[1:4 6:9]
        fprintf(tfile,['cmd_moveMotor_DurationInterval 1,' num2str(pid) ',' num2str(stage) ',' num2str(initSteps) ',' num2str(initOt) ',2.5\n']);
    end
    fprintf(tfile, [
        ['cmd_waitSeconds ' num2str(abs(initSteps)*2.5/1000+1) '\n']...
        'cmd_waitSeconds  2\n']);
    
    
end



fclose(tfile);
end

