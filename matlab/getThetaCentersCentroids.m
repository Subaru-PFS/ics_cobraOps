%% MATLAB CENTROID EXTRACTING SCRIPT


%% USER INPUTS AND DIRECTIONS

% CLEAR WORKSPACE IF USING THIS SCRIPT STANDALONE
clear all
clear data;
close all;
% Provide directory containing images
ImgDir = '~/temp/ten_images/';
if ~exist('ImgDir','var'), ImgDir = 'D:\PfsTests\15_04_15_15_13_00_msimMaps\Images/';end; %Must end with /
data.ImgDir = ImgDir;

dataName = ['s1cnrtrs_' char(regexp(pwd,'(\d+_)*','match'))];

% Specify a filter for image name (regexp)
if ~exist('ImgFilter','var'), ImgFilter = '*s1*.fits';end; % make sure to specify fits or ppm 
ImgFilter = 'theta*.fits';
data.ImgFilter = ImgFilter;

% Comment this out if you don't want to save the resulting structure
if ~exist('savedir','var'), savedir = ImgDir;end;

% Assign all centroids to one data structure or seperate structures for each
% image
% SET THIS FALSE FOR STROBED MODE IMGS
if ~exist('assgn2one','var'), assgn2one = false;end; % True for one structure. False for individual image structures.

% Which positioners are in the images? pId = (moduleId-1)*57+railPos 
% Example: Cobra 2.5 is pId 62
data.pId = [1,3,5,7,9,11,13,15,17,19,21,23,25,27,29,31,33,35,37,39,41,43,45,47,49,51,53]; % Temporary numbering for EM initial build

% Positioners of interest for center metrology. Which ones are moving? Use
% the array indecies from pId array.
data.activeCobras = [1,5,9,13,17,21,25,29,33,37,41,45,49,53]; % Same purpose as movingCobras but more intuitive to use pId

% Positioners which should be treated as fiducials. @Johannes change for
% debug
data.fixedCobras = [3,7,11,15,19,23,27,35,39,43,47,51,]; % Same purpose as Fiducials, but more intuitive to use pId

%Fiducal Numbers w/o fixed cobras.
data.fiducials = [101, 102, 103, 104, 105];

% Get the config xml for pulling centers
CobraCfg = loadCfgXml;

% Specify number of sigma above image median to threshold centroids at
% (leave at 5 if unsure)
NSIGMA = 3;

% Cobra Centers.  This determines the centroid assignments.  
% Use as accurate numbers as possible.  For fiducials, use the 
% actual location of the fiber, not the center of the fiducial cobra.  
% For moving positioners, use the center of S1.  These should be in 
% the same order as the pId array.
% The following arrays has centers for pid 1-9 in consecutive order
% PID 1,2,3...

%Module ID in your XML 
mID = 1;
c=1;
for pid=data.activeCobras
    data.centers(c) = getARMval(CobraCfg,pid, mID,'Global_base_pos_x')+getARMval(CobraCfg,pid, mID,'Global_base_pos_y')*i;
    data.Rptrl(c) = 60;
    c=c+1;
end
data.Rptrl(6) = 90;
for pid=data.fixedCobras
    data.centers(c) = getARMval(CobraCfg,pid, mID,'Global_base_pos_x')+getARMval(CobraCfg,pid, mID,'Global_base_pos_y')*i;
    data.Rptrl(c) = 7;
    c=c+1;
end
 
% Add the fiducial positions to centers array
data.centers = [data.centers.';...
72.7181901965125 + 362.439887112847i;...
461.498446356842 + 733.045426660864i;...
850.611264646654 + 1104.98112670259i;...
1238.52772825481 + 1476.34296141023i;...
1628.80387570267 + 1849.41980465647i];

% Control Region
data.Rptrl =  [data.Rptrl.';
    10;
    10;
    10;
    10;
    10];
   
%keyboard;
%% CENTROID PROCESSING (DONT CHANGE)

switch assgn2one
    
    case true
        % Get centroids from folder
        data = getCentroidsFromFolder(data);

        % Copy data structure with directory name
        assignin('base',dataName,data);

        if exist('savedir','var')
            save(strcat(savedir,dataName,'.mat'),dataName);
        end
        
    case false
        images = ls([data.ImgDir data.ImgFilter]);
        clear imgNames;
        
        for img=1:length(images(:,1))

            disp(strcat('Now processing ',images(img,:)));

            %Clear the temp data structure
            clear tdata;
            
            tdata = data;
            tdata.Name = char(images(img,:));
            tdata.imageName = images(img,:);

            % Get centroids and assign to correct Cobra or fiducial
            tdata = getCntrdsFromImage(tdata,NSIGMA);
        

            % Copy data structure with directory name
            imgName = char(regexp(tdata.Name,'.*ImageId_\d*','match'));
            assignin('base',imgName,tdata);
            imgNames{img} = imgName;
        end
        
        if exist('savedir','var')
            save(fullfile(savedir,strcat(dataName,'.mat')),imgNames{:});
        end
end



%% INTENSITY STABILITY PLOT
% figure
% CW = lines(8);
% for ii=1:8
%     fld = sprintf('pId%d',ii);
%     subplot(3,3,ii)
%     plot(data.(fld).pxSum,'color',CW(ii,:))
% end



