%% ON TIME SCRIPT GENERATION
% This makes the ontime tuning MSIM scripts using streaks
clear all

%% INPUTS
% Set exposure time
TexpTheta = 40;%sec
TexpPhi = 20;%sec
Tstrobe = 0.5;%sec   %Do not comment out, but not used when streaks are used
phiRefOt = 0.15;
thetaRefOt = .3;
pids = [1:9]; % Which pIds in MSIM will be commanded?
invertedThetas = [6]; % Which pids have theta motors which are inverted?
invertedPhis = [];

    % THETA ======================
    stage = 1;
    ontimes = 0.04:0.02:0.30; % Specify which ontimes to try
    
    cntrStrkOt = 0.2; % Ontime to use for finding centers in the beginning. Should be slow
    % enough for a good streak but strong enough to make it around. This is
    % a blind guess in the beginning and will have different results for
    % various positioners. A full circle is not necessary, just enough to
    % fit one and find the center.
    
    
    % Ontime to use for sending positioners back to home or antihome
    initOt = 0.3;
    
        % FORWARD SCRIPT
        filename = 'C:\Users\cmorantz\Dropbox\PFS_EM\SVN\MSIM\20140718_thetaFwdOntimeTuning.lst';
        thetaRefSteps = 7000; %steps to get to antihom for reference centroid
        phiRefSteps = 7000; % steps to get to antihome for reference centroid
        cntrStrkSteps = -7000; %steps to do center finding streak
        romSteps = 7000; % Steps to use for ontime tests
        initSteps = -8000; % Steps to command back home with
        makeOntimeScript(filename,pids,invertedThetas,invertedPhis,TexpTheta,Tstrobe,ontimes,thetaRefSteps,thetaRefOt,phiRefSteps,phiRefOt,stage,cntrStrkSteps,cntrStrkOt,romSteps,initSteps,initOt);
        
        % REVERSE SCRIPT
        filename = 'C:\Users\cmorantz\Dropbox\PFS_EM\SVN\MSIM\20140718_thetaRevOntimeTuning.lst';
        thetaRefSteps = -7000; %steps to go home on theta for reference centroid
        phiRefSteps = 7000; % steps to go home on phi for reference centroid
        cntrStrkSteps = 7000; %steps to do center finding streak
        romSteps = -7000;% Steps to use for ontime tests
        initSteps = 8000;% Stepst to command out to antihome to start tests from
        makeOntimeScript(filename,pids,invertedThetas,invertedPhis,TexpTheta,Tstrobe,ontimes,thetaRefSteps,thetaRefOt,phiRefSteps,phiRefOt,stage,cntrStrkSteps,cntrStrkOt,romSteps,initSteps,initOt);
        
    
    % PHI ======================
    stage=2;
    ontimes = 0.05:0.01:0.15;
    cntrStrkOt = 0.05;
    initOt = 0.15;
    
        % FORWARD SCRIPT
        filename = 'C:\Users\cmorantz\Dropbox\PFS_EM\SVN\MSIM\20140718_phiFwdOntimeTuning.lst';
        thetaRefSteps = 4000;
        phiRefSteps = 4000;
        cntrStrkSteps = -4000;
        romSteps = 3000;
        initSteps = -4000;
        makeOntimeScript(filename,pids,invertedThetas,invertedPhis,TexpPhi,Tstrobe,ontimes,thetaRefSteps,thetaRefOt,phiRefSteps,phiRefOt,stage,cntrStrkSteps,cntrStrkOt,romSteps,initSteps,initOt);
        
        % REVERSE SCRIPT
        filename = 'C:\Users\cmorantz\Dropbox\PFS_EM\SVN\MSIM\20140718_phiRevOntimeTuning.lst';
        thetaRefSteps = -4000;
        phiRefSteps = -4000;
        cntrStrkSteps = 4000;
        romSteps = -3000;
        initSteps = 4000;
        makeOntimeScript(filename,pids,invertedThetas,invertedPhis,TexpPhi,Tstrobe,ontimes,thetaRefSteps,thetaRefOt,phiRefSteps,phiRefOt,stage,cntrStrkSteps,cntrStrkOt,romSteps,initSteps,initOt);





