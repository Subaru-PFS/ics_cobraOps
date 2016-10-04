function output=getARMval(CobraConfig, pID, mID, parameter)
% get a numeric arm kinematics parameter
output = [];
ncobras = length(CobraConfig.ARM_DATA.ARM_DATA_CONTAINER);
%% generate arrays to help find container index for pid/mid
for jj=1:ncobras
    pid(jj) = str2num(CobraConfig.ARM_DATA.ARM_DATA_CONTAINER{jj}.DATA_HEADER.Positioner_Id.Text);
    mid(jj) = str2num(CobraConfig.ARM_DATA.ARM_DATA_CONTAINER{jj}.DATA_HEADER.Module_Id.Text);
end
    
if ~exist('pID','var')
    disp('Second parameter is a pID number:');
    disp(unique(pid));
    return
end
if ~exist('mID','var')
    disp('Third parameter is the module ID:');
    disp(unique(mid));
    return
end
%% get the container index for the pID,mID requested
kk = find(pID == pid & mID == mid,1);
if isempty(kk)
   return 
end
if ~exist('parameter','var')
    fprintf(1,'Valid field names for pid %d mid %d:\n',pID,mID);
    disp(fieldnames(CobraConfig.ARM_DATA.ARM_DATA_CONTAINER{kk}.KINEMATICS))
    return
end

%% is this conditional vestigial?
if(ischar(pID))
    pID = str2num(pID);
end 

%keyboard;
if(isfield(CobraConfig.ARM_DATA, 'ARM_DATA_CONTAINER')) % For backward compatibility reasons.
    if isfield(CobraConfig.ARM_DATA.ARM_DATA_CONTAINER{kk}.KINEMATICS, parameter)
        output = str2num(CobraConfig.ARM_DATA.ARM_DATA_CONTAINER{kk}.KINEMATICS.(parameter).Text);
        return
    elseif isfield(CobraConfig.ARM_DATA.ARM_DATA_CONTAINER{kk}, parameter)
        %% remove this condition when phiMin and phiMax are moved into ARM_DATA_CONTAINER
        output = str2num(CobraConfig.ARM_DATA.ARM_DATA_CONTAINER{kk}.(parameter).Text);
        return
    else
        warning(strcat(parameter, ' - Parameter name is not in configuration structure.'));
        %disp(fieldnames(CobraConfig.ARM_DATA.ARM_DATA_CONTAINER{kk}.KINEMATICS));
    end
else
    %% we should just dump this condition.
    for ii = 1:sum(regexpcmp(fieldnames(CobraConfig.ARM_DATA),'ARM_DATA_*')) % Count all struct that start with the name ARM_DATA_
        ARM_DATA_n = sprintf('ARM_DATA_%d', ii);
        if(pID ==  str2num(CobraConfig.ARM_DATA.(ARM_DATA_n).DATA_HEADER.Positioner_Id.Text))
            if isfield(CobraConfig.ARM_DATA.(ARM_DATA_n).KINEMATICS, parameter)
                output = str2num(CobraConfig.ARM_DATA.(ARM_DATA_n).KINEMATICS.(parameter).Text);
                return
            elseif isfield(CobraConfig.ARM_DATA.(ARM_DATA_n), parameter)
                output = str2num(CobraConfig.ARM_DATA.(ARM_DATA_n).(parameter).Text);
                return
            else
                warning(strcat(parameter, ' - Parameter name is not in configuration structure.'));
            end
        end
    end
end 



