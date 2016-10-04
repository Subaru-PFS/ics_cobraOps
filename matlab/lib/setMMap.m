function output=setMMap(CobraConfig, armID, speed, JointID, direction, data, data_regs)
% set FAST or SLOW motor map table from numeric array
% CobraConfig: xml configuration file read in by []
% armID      : numeric arm (positioner) ID
% speed      : 'FAST' or 'SLOW'
% JointID    : 1 or 2 (numeric), the joint ID number
% direction  : 'fwd' or 'rev'
% data       : the table in a numeric array
%            : 1x(n+2) for steps only (for backwards compatibility)
%            : 2x(n+2) for regions and steps (Johannes's input)
% data_regs  : optional stepsize array (additional functionality)
% returns an updated version of CobraConfig

output = CobraConfig;
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

TBL_SPEED   = [speed '_CALIBRATION_TABLE'];
TABLE_REGS  =  sprintf('Joint%d_%s_regions', JointID, direction);
TABLE_STEPS = sprintf('Joint%d_%s_stepsizes', JointID, direction);

if(isfield(CobraConfig.ARM_DATA, 'ARM_DATA_CONTAINER')) % For backward compatibility reasons. 
for jj = 1: length(CobraConfig.ARM_DATA.ARM_DATA_CONTAINER)
    if(armID ==  str2num(CobraConfig.ARM_DATA.ARM_DATA_CONTAINER{jj}.DATA_HEADER.Positioner_Id.Text))        
        if isfield(output.ARM_DATA.ARM_DATA_CONTAINER{jj}.(TBL_SPEED), TABLE_STEPS)
            if(size(data,2) >= 2)
                output.ARM_DATA.ARM_DATA_CONTAINER{jj}.(TBL_SPEED).(TABLE_REGS).Text  = NumArr2mMapStr(data(1,:));
                output.ARM_DATA.ARM_DATA_CONTAINER{jj}.(TBL_SPEED).(TABLE_STEPS).Text = NumArr2mMapStr(data(2,:));
            else
                if exist('data_regs','var')
                    output.ARM_DATA.ARM_DATA_CONTAINER{jj}.(TBL_SPEED).(TABLE_REGS).Text = NumArr2mMapStr(data_regs);
                end
                output.ARM_DATA.ARM_DATA_CONTAINER{jj}.(TBL_SPEED).(TABLE_STEPS).Text  = NumArr2mMapStr(data);
            end
        else
            warning('Parameter name does not exist in configuration structure');
        end
    end
end
else 
    ARM_DATA_n  = sprintf('ARM_DATA_%d', armID);
for i = 1:sum(regexpcmp(fieldnames(CobraConfig.ARM_DATA),'ARM_DATA_*')) % Count all struct that start with the name ARM_DATA_
    ARM_DATA_n = sprintf('ARM_DATA_%d', i);
    if(armID ==  str2num(CobraConfig.ARM_DATA.(ARM_DATA_n).DATA_HEADER.Positioner_Id.Text))
        
        if isfield(output.ARM_DATA.(ARM_DATA_n).(TBL_SPEED), TABLE_STEPS)
            if(size(data,2) >= 2)
                output.ARM_DATA.(ARM_DATA_n).(TBL_SPEED).(TABLE_REGS).Text  = NumArr2mMapStr(data(1,:));
                output.ARM_DATA.(ARM_DATA_n).(TBL_SPEED).(TABLE_STEPS).Text = NumArr2mMapStr(data(2,:));
            else
                if exist('data_regs','var')
                    output.ARM_DATA.(ARM_DATA_n).(TBL_SPEED).(TABLE_REGS).Text = NumArr2mMapStr(data_regs);
                end
                output.ARM_DATA.(ARM_DATA_n).(TBL_SPEED).(TABLE_STEPS).Text  = NumArr2mMapStr(data);
            end
        else
            warning('Parameter name does not exist in configuration structure');
        end
    end
end
end
 

