% TL = dlmread('C:\Users\cmorantz\Dropbox\PFSMatlab\Results\FiducialMeasurements\02_20_14_12_55_23_CentroidRtblty\Centroid_log.txt',',');
function output = processCentroidLog(logFileName)


logFile = fopen(logFileName);

% Read first log file line
thisLine = fgetl(logFile);

% Determine how many cobras there are and initialize arrays
pid = str2num(char(regexp(thisLine,'(?<=Pid=)[^,]*','match')));
mid = str2num(char(regexp(thisLine,'(?<=Mid=)[^\.]*','match')));
nPid = length(pid);
%for ii=1:nPid

for ii = 1:length(pid)
    % Label positioners sequenctially starting at pId1
    fldId = strcat('pId',num2str(pid(ii)));
    % OR 
    % Label positioners with their PFI index
    %       fldId = strcat('pId',num2str((mid(ii)-1)*57+pid(ii)));
    
    % Initialize array
    data.(fldId).CLpos = [];
    data.(fldId).pid = pid(ii);
end

% Determine how many fiducials there are and initialize arrays
fid = str2num(char(regexp(thisLine,'(?<=Fid=)[^,]*','match')));
nFid = length(fid);
data.nFid = nFid;
for ii=1:nFid
    % Start fiducial pIds counting where positioner pIds left off
    fldId = strcat('fId',num2str(ii));
    data.(fldId).CLpos = [];
end

while thisLine ~= -1
    
    % Pull out fiducial and positioner IDs from line
    fid = str2num(char(regexp(thisLine,'(?<=Fid=)[^,]*','match')));
    mid = str2num(char(regexp(thisLine,'(?<=Mid=)[^\.]*','match')));
    pid = str2num(char(regexp(thisLine,'(?<=Pid=)[^,]*','match')));
    
    % Check that line contains the correct number of entries
    if (length(fid)~=nFid) | (length(pid)~=nPid)
        disp('throwing out line')
        disp(thisLine)
        thisLine = fgetl(logFile);
        continue;
    end
    
    % Pull out centroids from line
    X = str2num(char(regexp(thisLine,'(?<=X=)[^,]*','match')));
    Y = str2num(char(regexp(thisLine,'(?<=Y=)[^ ]*','match')));

    % Save fiducial positions to CLpos array
    for ii=1:nFid
        fldId = strcat('fId',num2str((ii)));
        data.(fldId).CLpos = [data.(fldId).CLpos, X((ii)) + Y((ii))*i];
    end
    
    % Save cobra positions to CLpos array
    for ii=1:nPid
        fldId = strcat('pId',num2str(pid(ii)));
        data.(fldId).CLpos = [data.(fldId).CLpos, X(nFid+(ii))+Y(nFid+(ii))*i];
    end
    
    thisLine = fgetl(logFile);

end

output = data;

return

