% Get ROM of Cobras
clear all
close all

%% INPUTS

% Provide mat file containing the structures from processing the ROM and
% orientation images.
[romImgMat romImgMatPath] = uigetfile('..\..\TEST_RESULTS\Metrology\*.mat');

% Provide the mat file for array of reference fiducial positions
[refposFile refposFilePath] = uigetfile('*.mat','Choose reference fiducial position file');

% Provide config file with link lengths
% configStructFile = 'CobraConfig_061014_5_inv.mat';
% configFileXml = 'EMcfg_StrkAdj1_061714.xml';
CobraConfig = loadCfgXml;

%% Calculate ROM and orientation
% Load img structures
load(fullfile(romImgMatPath,romImgMat));
load(fullfile(refposFilePath,refposFile));


S = whos;

% THETA IN
A = regexp(cellstr([S.name]),'thetaInPhiOut_ImageId_\d*','match');
thetaIn = eval(A{1}{1});
existingFigure('Orientation Reference img')
img = fitsread('..\..\dotsBefore_ImageId_3_Target_0_Iteration_-1_loopId_-1_64.fits');
imgflipped = flipdim(img,1);
imagesc(imgflipped)
axis equal
hold on
thetaIn = applyHorns(thetaIn,refpos);
thetaIn = imgJointAngles(thetaIn,CobraConfig,true,true);

% THETA OUT
A = regexp(cellstr([S.name]),'thetaPhiOut_ImageId_\d*','match');
thetaOut = eval(A{1}{1});
thetaOut = applyHorns(thetaOut,refpos);
thetaOut = imgJointAngles(thetaOut,CobraConfig,true,false);

activeCobras = thetaIn.activeCobras;

for ii = activeCobras
    fldID = sprintf('pId%d',ii);
    cnfgID = sprintf('ARM_DATA_%d',ii);
    
    CobraConfig.ARM_DATA.(cnfgID).KINEMATICS.Global_base_ori_z.Text = num2str(thetaIn.(fldID).J1*180/pi);
    
    CobraConfig.ARM_DATA.(cnfgID).KINEMATICS.Joint1_CCW_limit_angle.Text = '0';
    CobraConfig.ARM_DATA.(cnfgID).KINEMATICS.Joint2_CCW_limit_angle.Text = '0';

    % Grab centroids from theta in/out images
    thetaInJ1 = thetaIn.(fldID).J1;
    thetaOutJ1 = thetaOut.(fldID).J1;
   
    % Check that thetaOut is not overtaking thetaIn and crossing +X axis
    % Add 2*pi if it is
    if (thetaOutJ1 - thetaInJ1) < -3
        thetaOutJ1 = thetaOutJ1 + 2*pi;
        disp(strcat('added 2pi to thetaOut on ',fldID))
    end

    % Calculate theta rom assuming that thetaOut overtakes thetaIn
    rom = 2*pi + (thetaOutJ1 - thetaInJ1);
    
%     CobraConfig.ARM_DATA.(cnfgID).KINEMATICS.Joint1_CW_limit_angle.Text = num2str(rom*180/pi);
    CobraConfig.ARM_DATA.(cnfgID).KINEMATICS.Joint1_CW_limit_angle.Text = '400';
    CobraConfig.ARM_DATA.(cnfgID).KINEMATICS.Joint2_CW_limit_angle.Text = num2str(min(thetaIn.(fldID).J2,thetaIn.(fldID).J2)*180/pi);
    
    
end


