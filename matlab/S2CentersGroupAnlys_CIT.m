    % S1CentersGroupAnlys

clear all
close all

c=1;
imgsFiles = {};
dataTitle = {};
mID = 1; % Module ID in XML;
CobraConfig = loadCfgXml;

PxSc = 90.5;

metDir = '.';
metDirFiles = dir2cell(metDir);

phifltr = 's2cntrs_.*';

phiFiles = regexp(metDirFiles,phifltr,'match');
phiFiles = phiFiles(~cellfun('isempty',phiFiles));
if length(phiFiles)==1
    imgsFiles = {fullfile(metDir,char(phiFiles{1}))};
else
    imgsFiles = fullfile(metDir,[phiFiles{:}]);
end
dataTitle = [dataTitle phiFiles];
 
 
%% EXECUTION

Link1s = [];
Link2s = [];
dates = [];
L1Means = [];
L2Means = [];

for ii=1:length(imgsFiles(:,1))
    
  
     data{ii} = S2CentersAnlys_CIT(imgsFiles{ii},dataTitle{ii},CobraConfig, mID);
 
    
    Link1s = vertcat(Link1s, data{ii}.Link1s);
    Link2s = vertcat(Link2s, data{ii}.Link2s);
 
    L1Means = vertcat(L1Means, mean(data{ii}.Link1s(data{ii}.Link1s>20)));
    L2Means = vertcat(L2Means, mean(data{ii}.Link1s(data{ii}.Link2s<60 & data{ii}.Link1s>20))); % Filter 
    

end

keyboard;

hf3 = existingFigure('Link Lengths');
numPos = length(Link1s(1,:));
for ii=1:numPos
    TLA = Link1s(:,ii);% [Temporary Link Array]
    TLA(TLA==0) = []; % Deleting zeros
    if(~isempty(TLA))
    tlMean = mean(TLA); % Take mean of links. Filtering out those greater than 60px
    CobraConfig = setARMval(CobraConfig,data{1}.activeCobras(ii), mID,'Link1_Link_Length',tlMean);
    tlStd = std(TLA); % Take standard deviation of links. Filtering out those greater than 60px
    figure(hf3)
    subplot(numPos,2,2*(ii-1)+1)
 legend(sprintf('Overall Mean=%4.2fum \\sigma=%4.2fum',tlMean*PxSc,tlStd*PxSc));
    end
end

numPos = length(Link2s(1,:));
for ii=1:numPos
    TLA = Link2s(:,ii);% [Temporary Link Array]
    TLA(TLA==0) = []; 
      if(~isempty(TLA))
    tlMean = mean(TLA); % Take mean of links. Filtering out those greater than 60px
    CobraConfig = setARMval(CobraConfig,data{1}.activeCobras(ii), mID,'Link2_Link_Length',tlMean);
    tlStd = std(TLA); % Take standard deviation of links. Filtering out those greater than 60px
    figure(hf3)
    subplot(numPos,2,2*ii)
    legend(sprintf('Overall Mean=%4.2fum \\sigma=%4.2fum',tlMean*PxSc,tlStd*PxSc));
      end
end


[xmlfile, xmlfilepath] = uiputfile('*.xml','Save new CobraConfig XML file with centers');
cobraCfg2xml(CobraConfig,fullfile(xmlfilepath,xmlfile));  
