function CobraConfig=setARMval(CobraConfig, pID, mID, parameter, value)
% set a numeric or text arm kinematics parameter.
% parameter must ALREADY exist in CobraConfig
% CobraConfig = config structure
% armID = integer ID number of positioner
% parameter = field name to be updated (string)
% value = value of paramter

%% check inputs
if ~exist('pID','var')
    disp('Second parameter is an pID number');
     return
end
if ~exist('mID','var')
    disp('Third parameter is the module ID');
     return
end
if ~exist('parameter','var')
    disp('CAUTION THIS METHOD HAS CHANGED. From using ARM_ID to pID and module ID. Fourth parameter is a field name from the KINEMATICS structure:');
    return
end

%% set defaults
if ~exist('value','var')
        disp('CAUTION THIS METHOD HAS CHANGED. From using ARM_ID to pID and module ID. Fifth parameter is the value you want to set');
    disp('No value for parameter set - using empty string');
    value = '';
end

%% convert numeric value to string
if isnumeric(value)
    inputvalue = num2str(value);
elseif ischar(value)
    inputvalue = value;
else
    warning('class of value should be numeric or char');
    return;
end
 

if(isfield(CobraConfig.ARM_DATA, 'ARM_DATA_CONTAINER')) % else case for backward compatibility reasons. 
    for jj = 1: length(CobraConfig.ARM_DATA.ARM_DATA_CONTAINER)
    if(pID ==  str2num(CobraConfig.ARM_DATA.ARM_DATA_CONTAINER{jj}.DATA_HEADER.Positioner_Id.Text))       
       % if isfield(CobraConfig.ARM_DATA.ARM_DATA_CONTAINER{jj}.KINEMATICS, parameter) 
            CobraConfig.ARM_DATA.ARM_DATA_CONTAINER{jj}.KINEMATICS.(parameter).Text = ...
                inputvalue
            
            return;
      %  elseif isfield(CobraConfig.ARM_DATA.ARM_DATA_CONTAINER{jj}, parameter)
      %       CobraConfig.ARM_DATA.ARM_DATA_CONTAINER{jj}.(parameter).Text = ...
      %          inputvalue;
      %  elseif strcmp(parameter, 'phiMax')| strcmp(parameter, 'phiMin')
      %       CobraConfig.ARM_DATA.ARM_DATA_CONTAINER{jj}.(parameter).Text = ...
      %          inputvalue;
      %  else
      %      warning(strcat(parameter, ' Parameter name is not in configuration structure. Does not get stored in XML!'));
      %  end
    end
end
else
for i = 1:sum(regexpcmp(fieldnames(CobraConfig.ARM_DATA),'ARM_DATA_*')) % Count all struct that start with the name ARM_DATA_
    ARM_DATA_n = sprintf('ARM_DATA_%d', i);
    if(pID ==  str2num(CobraConfig.ARM_DATA.(ARM_DATA_n).DATA_HEADER.Positioner_Id.Text))
        if isfield(CobraConfig.ARM_DATA.(ARM_DATA_n).KINEMATICS, parameter) 
            CobraConfig.ARM_DATA.(sprintf('ARM_DATA_%d', ...
                i)).KINEMATICS.(parameter).Text = ...
                inputvalue;
            return;
        elseif isfield(CobraConfig.ARM_DATA.(ARM_DATA_n), parameter)
             CobraConfig.ARM_DATA.(sprintf('ARM_DATA_%d', ...
                i)).(parameter).Text = ...
                inputvalue;
        else
            warning('Parameter name is not in configuration structure');
        end
    end
end
end
