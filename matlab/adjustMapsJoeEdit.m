%% DAMPING MODIFICATION USING TARGET CONVERGENCE DATA
clear all
close all

% Load mat files in current directory
mIdMats = loadmats('*mId*pId*.mat');
eval(unpackstruct(mIdMats));

% Load cobra config xml
[cfgFile,cfgPath] = uigetfile('*.xml','Select the original XML config file used for this test');
CobraConfig = xml2struct(fullfile(cfgPath,cfgFile));

% Specify which maps will be updated (rev, fwd or both)
revUpdate = true;
fwdUpdate = true;

% Specify number of target iterations to use (should make this greater than what
% is allowed by MSIM script by default)
maxSteps = 16;

% ---------------------------END OF INPUT SECTION-------------------------

% Store current CobraConfig struct
CobraConfigOrig = CobraConfig;
clear CobraConfig;
% Run modify damping script
CobraConfig = modifyDampingJoeEdit(CobraConfigOrig, '.', maxSteps, revUpdate, fwdUpdate);
% Save new XML with adjusted maps
% [xmlfile, xmlfilepath] = uiputfile('*.xml','Save new CobraConfig XML file with adjusted motor maps');
% cobraCfg2xml(CobraConfig,fullfile(xmlfilepath,xmlfile));