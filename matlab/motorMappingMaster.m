%% MOTOR MAPPING MASTER
% This script should be run from    Dropbox\PFS_EM\SVN\MATLAB

%% OPTIONAL CLEAR AND CLOSE
clear all
close all

%% STORE STREAK MAPS TO COBRA CONFIG STRUCTURE
if false % RUN ME OR NOT?
% Give initial CobraConfig structure to which the motor maps will be saved
CobraConfigOrig = 'CobraConfig_061014_5_inv.mat';
% Give the locations of motor map workspaces output by the streak motor map
% autogeneration script.
thetaMMapWorkspace = '..\..\TEST_RESULTS\StreakResults\workspacestage1.mat';
phiMMapWorkspace = '..\..\TEST_RESULTS\StreakResults\workspacestage2.mat';
% Number of theta bins in motor map
nTheta = 100;
% Number of phi bins in motor map
nPhi = 45;
% pIdMapping, see streakMaps2cfgStrct source for description. Comment out
% if not necessary
pIdMapping = [1,2,3,4,7,8,9];
% ---------------------------END OF INPUT SECTION-------------------------
% Store motor maps to structure
load(CobraConfigOrig);
CobraConfigOrig = CobraConfig;
clear CobraConfig;
CobraConfig = streakMaps2cfgStrct(CobraConfigOrig, thetaMMapWorkspace, phiMMapWorkspace, nTheta, nPhi, pIdMapping);
end


%% DAMPING MODIFICATION USING TARGET CONVERGENCE DATA
if true % RUN ME OR NOT?
% Provide initial Cobra config matlab structure with motor maps to be
% modified
cnfgFile = 'CobraConfig_061614_1.mat';
% Provide target convergence data directory for initial maps
dataDir = '..\..\TEST_RESULTS\TargetConvergence\06_17_14_11_21_32_TargetRun';
% Specify which maps will be updated (rev, fwd or both)
revUpdate = false;
fwdUpdate = true;
% Specify number of target iterations to use (should make this greater than what
% is allowed by MSIM script by default)
maxSteps = 16;
% ---------------------------END OF INPUT SECTION-------------------------
% Load config file
clear CobraConfig;
load(cnfgFile);
CobraConfigOrig = CobraConfig;
clear CobraConfig;
% Run modify damping script
CobraConfig = modifyDamping(CobraConfigOrig, dataDir, maxSteps, revUpdate, fwdUpdate);
end
