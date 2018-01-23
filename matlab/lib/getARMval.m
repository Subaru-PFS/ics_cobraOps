function output=getARMval(CobraConfig, pID, mID, parameter)
% get a numeric arm kinematics parameter(s)
%
% Usage: output=getARMval(CobraConfig, pID, mID, parameter)
% 
% CobraConfig: data structure ready in by loadCfgXml.m
% pID: positioner ID - may be an array
% mID: module ID.  Always 1 for now, but in the future may be
% something else
% parameter: char array (old string)  argument name of parameter
%
% missing parameters lead to a listing of acceptable values and
% null return value.

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
%% get the container indices for the [pID],mID requested
mlocations = find(mID == mid);
%kk = find(pID == pid & mID == mid,1);
[junk_logical plocations] = ismember(pID,pid(mlocations)); % plocations are only the first location
                                                           % where it matches, so we have to
                                                           % filter on mids that match the request.
containers = mlocations(plocations); % yes this is convoluted, but it allows for multiple values of mid.
if isempty(containers)
   return 
end
if ~exist('parameter','var')
    fprintf(1,'Valid field names for pid %d mid %d:\n',pID,mID);
    disp(fieldnames(CobraConfig.ARM_DATA.ARM_DATA_CONTAINER{containers(1)}.KINEMATICS))
    return
end

%% is this conditional vestigial?
if(ischar(pID))
    pID = str2num(pID);
end 

if(isfield(CobraConfig.ARM_DATA, 'ARM_DATA_CONTAINER')) % For backward compatibility reasons.
    for jj = 1:length(containers)
        kk = containers(jj);
        if isfield(CobraConfig.ARM_DATA.ARM_DATA_CONTAINER{kk}.KINEMATICS, parameter)
            output(jj) = str2num(CobraConfig.ARM_DATA.ARM_DATA_CONTAINER{kk}.KINEMATICS.(parameter).Text);
        elseif isfield(CobraConfig.ARM_DATA.ARM_DATA_CONTAINER{kk}, parameter)
            %% remove this condition when phiMin and phiMax are moved into ARM_DATA_CONTAINER
            output(jj) = str2num(CobraConfig.ARM_DATA.ARM_DATA_CONTAINER{kk}.(parameter).Text);
        else
            warning(strcat(parameter, ' - Parameter name is not in configuration structure.'));
            %disp(fieldnames(CobraConfig.ARM_DATA.ARM_DATA_CONTAINER{kk}.KINEMATICS));
        end
    end
else
    disp('Configuration structure appears malformed');
end


