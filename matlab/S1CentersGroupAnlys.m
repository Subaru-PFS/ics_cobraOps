% S1CentersGroupAnlys

clear all
close all

c=1;
imgsFiles = {};
dataTitle = {};
configFilePath = 'EMconfig_06-05-14.mat';
configStructPath = 'CobraConfig_060914_4.mat';
configSavePath = 'CobraConfig_060914_5.mat';

%3-26-14 S1 Centers
% imgsFiles{c} = '..\..\TEST_RESULTS\Metrology\centerms1_033114_1.mat';
% dataTitle{c} = '3/26/14-1';
% c=c+1;

metDir = '..\..\TEST_RESULTS\Metrology\';
metDirFiles = dir2cell(metDir);
thetaFileInd = strmatch('thetaCenters',metDirFiles);
thetaFiles = metDirFiles(thetaFileInd);
if length(thetaFiles)==1
    imgsFiles = [imgsFiles fullfile(metDir,thetaFiles{1})];
else
    imgsFiles = [imgsFiles fullfile(metDir,thetaFiles)];
end
dataTitle = [dataTitle thetaFiles];





load(configFilePath);
for ii=1:length(imgsFiles)
    data(ii) = S1CentersAnlys(imgsFiles{ii},dataTitle{ii},configFilePath);
end




% Calculate mean centers and save to xml structure
centersMean = mean(reshape([data.cBar].',[],length(data)).');

% Calculate overall variance
% http://en.wikipedia.org/wiki/Pooled_variance
sNmrtr = 0;
sDnmntr = 0;
for ii=1:length(data)
    tn = length(data{ii}.dnames);
    sNmrtr = sNmrtr + (tn-1)*(data{ii}.cSTD/sqrt(2)).^2;
    sDnmntr = sDnmntr + (tn-1);
end
centersStd = sqrt(sNmrtr/sDnmntr);


% if ~exist('CobraConfig','var')
%     CobraConfig = struct;
% end

load(configStructPath)

for ii=1:length(data(1).activeCobras)
    fldID = strcat('ARM_DATA_',num2str(data(1).activeCobras(ii)));
    CobraConfig.ARM_DATA.(fldID).KINEMATICS.Global_base_pos_x = real(centersMean(ii));
    
    % SUBTRACT Y FROM 2048 IF IMAGE FLIPPED FROM MSIM
%     CobraConfig.ARM_DATA.(fldID).KINEMATICS.Global_base_pos_y = 2048-imag(centersMean(ii));
    CobraConfig.ARM_DATA.(fldID).KINEMATICS.Global_base_pos_y = imag(centersMean(ii));

    CobraConfig.ARM_DATA.(fldID).KINEMATICS.Global_base_pos_z = 0;
end

save(configSavePath,'CobraConfig')

% cobraCfg2xml(CobraConfig, configSavePath);
    
