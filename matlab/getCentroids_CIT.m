%% MATLAB CENTROID EXTRACTING SCRIPT


%% USER INPUTS AND DIRECTIONS

% CLEAR WORKSPACE IF USING THIS SCRIPT STANDALONE
clear all
clear data;

% Data name
if ~exist('dataName','var'), dataName = 's1centers_062614';end;

% Provide directory containing images
if ~exist('ImgDir','var'), ImgDir = 'D:\PfsTests\06_26_14_10_19_55_centers\Images/';end; %Must end with /
data.ImgDir = ImgDir;

% Specify a filter for image name (regexp)
if ~exist('ImgFilter','var'), ImgFilter = '*s2*.fits';end; % make sure to specify fits or ppm 
data.ImgFilter = ImgFilter;

% Comment this out if you don't want to save the resulting structure
if ~exist('savedir','var'), savedir = ImgDir;end;

% Assign all centroids to one data structure or seperate structures for each
% image
% SET THIS FALSE FOR STROBED MODE IMGS
if ~exist('assgn2one','var'), assgn2one = false;end; % True for one structure. False for individual image structures.

% Which positioners are in the images? pId = (moduleId-1)*57+railPos 
% Example: Cobra 2.5 is pId 62
data.pId = [1,2,3,4,6,7,8,9,10,11,12]; % Temporary numbering for EM initial build

% Positioners of interest for center metrology. Which ones are moving? Use
% the array indecies from pId array.
data.activeCobras = [1,2,3,4,6,7,8,9]; % Same purpose as movingCobras but more intuitive to use pId

% Positioners which should be treated as fiducials. 
data.fixedCobras = [10,11,12]; % Same purpose as Fiducials, but more intuitive to use pId

% Get the config xml for pulling centers
CobraCfg = loadCfgXml;

% Cobra Centers.  This determines the centroid assignments.  
% Use as accurate numbers as possible.  For fiducials, use the 
% actual location of the fiber, not the center of the fiducial cobra.  
% For moving positioners, use the center of S1.  These should be in 
% the same order as the pId array.
% The following arrays has centers for pid 1-9 in consecutive order
% PID 1,2,3...

c=1;
for pid=data.activeCobras
    data.centers(c) = getARMval(CobraCfg,pid,'Global_base_pos_x')+getARMval(CobraCfg,pid,'Global_base_pos_y')*i;
    data.Rptrl(c) = 80;
    c=c+1;
end

% Add the fiducial positions to centers array
data.centers = [data.centers.';
    10^3*(1.6441081542969 + (2.048-0.1428441162109)*i);
    10^3*(1.2572365722656 + (2.048-0.5177436523438)*i);
    10^3*(0.0947736587524 + (2.048-1.6428116455078)*i)];

data.Rptrl =  [data.Rptrl.';
    10;
    10;
    10];
   

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
            tdata = getCntrdsFromImage(tdata);
        

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



