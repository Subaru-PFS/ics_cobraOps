% S1CentersGroupAnlys

clear all
close all

c=1;
imgsFiles = {};
dataTitle = {};

CobraConfig = loadCfgXml;

metDir = '.';
metDirFiles = dir2cell(metDir);
moduleID = 1;
thetafltr = 's1cnrtrs_.*';
 
thetaFiles = regexp(metDirFiles,thetafltr,'match');
thetaFiles = thetaFiles(~cellfun('isempty',thetaFiles));
if length(thetaFiles)==1
    imgsFiles = {fullfile(metDir,char(thetaFiles{1}))};
else
    imgsFiles = fullfile(metDir,[thetaFiles{:}]);
end
dataTitle = [dataTitle thetaFiles];


% Which positioners are in the images? pId = (moduleId-1)*57+railPos
% Example: Cobra 2.5 is pId 62
%data.pId = [1,3,5,7,9,11,13,15,17,19,21,23,25,27,29,31,33,35,37,39,41,43,45,47,49,51,53]; % Temporary numbering for EM initial build

% Positioners of interest for center metrology. Which ones are moving? Use
% the array indecies from pId array.
%data.activeCobras = [1,5,9,13,17,21,25,29,33,37,41,45,49,53]; % Same purpose as movingCobras but more intuitive to use pId

% Positioners which should be treated as fiducials.
%data.fixedCobras = [3,7,11,15,19,23,27,31,35,39,43,47,51,]; % Same purpose as Fiducials, but more intuitive to use pId
 
for ii=1:length(imgsFiles)
    data(ii) = S1CentersAnlys_CIT(imgsFiles{ii},dataTitle{ii},CobraConfig);
end
% Calculate mean centers and save to xml structure
if length(data)>1
    centersMean = mean(reshape([data.cBar].',[],length(data)).');
else
    centersMean = data.cBar;
end

keyboard;
% Calculate overall variance
% http://en.wikipedia.org/wiki/Pooled_variance
sNmrtr = 0;
sDnmntr = 0;
for ii=1:length(data)
    tn = length(data(ii).dnames);
    sNmrtr = sNmrtr + (tn-1)*(data(ii).cSTD/sqrt(2)).^2;
    sDnmntr = sDnmntr + (tn-1);
end
centersStd = sqrt(sNmrtr/sDnmntr);

for ii=1:length(data.activeCobras) % @Johannes
    pID =data.activeCobras(ii);
    if( length(centersMean)>=ii)
        if(real(centersMean(ii)) ~= 0) % we dont want to set it if for some reason its zero.
            CobraConfig = setARMval(CobraConfig,pID, moduleID, 'Global_base_pos_x',real(centersMean(ii)));
            CobraConfig = setARMval(CobraConfig,pID, moduleID, 'Global_base_pos_y',imag(centersMean(ii)));
        end
    end
end

[xmlfile, xmlfilepath] = uiputfile('*.xml','Save new CobraConfig XML file with centers');
cobraCfg2xml(CobraConfig,fullfile(xmlfilepath,xmlfile));

% Save figures in current directory
csfh = existingFigure('Center Spreads');
hgsave(csfh,'centerSpreads.fig');
csfh = existingFigure('Fiducial Spreads');
hgsave(csfh,'fidSpreads.fig');
