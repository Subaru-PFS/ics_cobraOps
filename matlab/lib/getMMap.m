function [outputY outputX]=getMMap(CobraConfig, armID, speed, JointID, direction)
% [Y X] = getMMap(CobraConfig, armID, speed, JointID, direction)
% get FAST or SLOW motor map table
% CobraConfig: xml configuration file read in by []
% armID      : numeric arm (positioner) ID
% speed      : 'FAST' or 'SLOW'
% JointID    : 1 or 2 (numeric), the joint ID number
% direction  : 'fwd' or 'rev'
% returns the table in a numeric array
%
% 9/16/14: PHM modified this to pull angles out while remaining
% backwards compatible with existing usage.
%
%output = [];
die = false;
if ~exist('armID','var') || ~isnumeric(armID)
    disp('Second parameter is an ARM_DATA_ ID number [1-9]');
    fieldnames(CobraConfig.ARM_DATA);                             die = true;
end
if ~exist('speed','var') || ~ischar(speed)
    disp('Third parameter is ''FAST'' or ''SLOW''');              die = true;
end
if ~exist('JointID','var') || ~isnumeric(JointID)
    disp('Fourth parameter is the Joint ID [1 or 2]');            die = true;
end
if ~exist('direction','var') || ~ischar(direction);
    disp('Fifth parameter is the direction ''fwd'' or ''rev''');  die = true;
end

if die, return; end;

%% check inputs for validity

if(ischar(armID))
    armID = str2num(armID);
end
%ARM_DATA_n = sprintf('ARM_DATA_%d', armID);
TBL_SPEED  = [speed '_CALIBRATION_TABLE'];
TABLE_DAT  = sprintf('Joint%d_%s_stepsizes', JointID, direction);
TABLE_ANG  = sprintf('Joint%d_%s_regions', JointID, direction);
if(isfield(CobraConfig.ARM_DATA, 'ARM_DATA_CONTAINER')) % For backward compatibility reasons.
for jj = 1: length(CobraConfig.ARM_DATA.ARM_DATA_CONTAINER)
    if(armID ==  str2num(CobraConfig.ARM_DATA.ARM_DATA_CONTAINER{jj}.DATA_HEADER.Positioner_Id.Text))
        if isfield(CobraConfig.ARM_DATA.ARM_DATA_CONTAINER{jj}.(TBL_SPEED), TABLE_DAT)
            outputY = mMapStr2NumArr(CobraConfig.ARM_DATA.ARM_DATA_CONTAINER{jj}.(TBL_SPEED).(TABLE_DAT).Text);
            outputX = mMapStr2NumArr(CobraConfig.ARM_DATA.ARM_DATA_CONTAINER{jj}.(TBL_SPEED).(TABLE_ANG).Text);
            outputX(3:end) = outputX(3:end)*pi/180;
            break;
        else
            warning('Parameter name does not exist in configuration structure');
        end
    end
end
else  % Backward compatibility when the xml was still havein ARM_DATA_1 etc.
for i = 1:sum(regexpcmp(fieldnames(CobraConfig.ARM_DATA),'ARM_DATA_*')) % Count all struct that start with the name ARM_DATA_
    ARM_DATA_n = sprintf('ARM_DATA_%d', i);
    if(armID ==  str2num(CobraConfig.ARM_DATA.(ARM_DATA_n).DATA_HEADER.Positioner_Id.Text))
        if isfield(CobraConfig.ARM_DATA.(ARM_DATA_n).(TBL_SPEED), TABLE_DAT)
            outputY = mMapStr2NumArr(CobraConfig.ARM_DATA.(ARM_DATA_n).(TBL_SPEED).(TABLE_DAT).Text);
            outputX = mMapStr2NumArr(CobraConfig.ARM_DATA.(ARM_DATA_n).(TBL_SPEED).(TABLE_ANG).Text);
            outputX(3:end) = outputX(3:end)*pi/180;
            break;
        else
            warning('Parameter name does not exist in configuration structure');
        end
    end
end
end