    % S1CentersGroupAnlys

clear all
close all

c=1;
imgsFiles = {};
dataTitle = {};
configFilePath = 'CobraConfig_060914_2.mat';
load(configFilePath);

% %3-28-14 S2 Centers
% imgsFiles{c} = '..\..\TEST_RESULTS\Metrology\centerms2_033114_1.mat';
% dataTitle{c} = '3/28/14-1';
% c=c+1;
% 
% metDir = '..\..\TEST_RESULTS\Metrology\';
% metDirFiles = dir2cell(metDir);
% phiFileInd = strmatch('phiMtrlgy',metDirFiles);
% phiFiles = {metDirFiles{phiFileInd}};
% imgsFiles = [imgsFiles fullfile(metDir,phiFiles)];
% dataTitle = [dataTitle phiFiles];

metDir = '..\..\TEST_RESULTS\Metrology\';
metDirFiles = dir2cell(metDir);
phiFileInd = strmatch('phiMtrlgy',metDirFiles);
phiFiles = metDirFiles(phiFileInd);
if length(phiFiles)==1
    imgsFiles = [imgsFiles fullfile(metDir,phiFiles{1})];
else
    imgsFiles = [imgsFiles fullfile(metDir,phiFiles)];
end
dataTitle = [dataTitle phiFiles];

%% EXECUTION

Link1s = [];
Link2s = [];
dates = [];
L1Means = [];
L2Means = [];


fh1 = figure('name','Link1s vs DOY');
fh2 = figure('name','Link2s vs DOY');
for ii=1:length(imgsFiles)
    data{ii} = S2CentersAnlys(imgsFiles{ii},dataTitle{ii},CobraConfig);
    Link1s = vertcat(Link1s, data{ii}.Link1s);
    Link2s = vertcat(Link2s, data{ii}.Link2s);
    
%     %Filter outliers
%     NN = length(data{ii}.Link1s(1,:));
%     for ll=1:NN
%         TD = data{ii}.Link1s(NN,:);
%         data{ii}.Link1s(NN,:) = TD(TD<60 & TD>20);
%     end
    
    L1Means = vertcat(L1Means, mean(data{ii}.Link1s(data{ii}.Link1s>20)));
    L2Means = vertcat(L2Means, mean(data{ii}.Link1s(data{ii}.Link2s<60 & data{ii}.Link1s>20)));
    
    tdate = datenum(regexp(imgsFiles{ii},'\d+_\d+_\d+','match'),'mm_dd_yy');
	dates = [dates tdate];
    
    
    figure(fh1)
    for jj = 1:length(data{ii}.Link1s(1,:))
        subplot(3,3,jj)
        plot(tdate*ones(length(data{ii}.Link1s),1),data{ii}.Link1s,'ro')
        hold on
    end
    
    figure(fh2)
    for jj = 1:length(data{ii}.Link2s(1,:))
        subplot(3,3,jj)
        plot(tdate*ones(length(data{ii}.Link2s),1),data{ii}.Link2s,'bo')
        hold on
    end
    
end

% figure(fh1)
% for jj = 1:
%     subplot(3,3,jj)
%     sortNplot(dates,L1Means(:,jj))
%     hold on
% end


hf3 = existingFigure('Link Lengths');
for ii=1:length(Link1s(1,:))
    TLA = Link1s(:,ii);% Three letter acronym meaning temporary link array
    tlMean = mean(TLA); % Take mean of links. Filtering out those greater than 60px
    
    fldID = strcat('ARM_DATA_',num2str(data{1}.activeCobras(ii)));
    CobraConfig.ARM_DATA.(fldID).KINEMATICS.Link1_Link_Length = tlMean;
    tlStd = std(TLA); % Take standard deviation of links. Filtering out those greater than 60px
    figure(hf3)
    subplot(9,2,2*(ii-1)+1)
    rlh = refline(0,tlMean,0,'k--');
%     legend(rlh,strcat('Overall Mean=',num2str(tlMean),' (\sigma=',num2str(tlStd),')'));
    legend(rlh,sprintf('Overall Mean=%4.2f \\sigma=%4.2f',tlMean,,tlStd));
end

for ii=1:length(Link2s(1,:))
    TLA = Link2s(:,ii);% Three letter acronym meaning temporary link array
    tlMean = mean(TLA); % Take mean of links. Filtering out those greater than 60px
    
    fldID = strcat('ARM_DATA_',num2str(data{1}.activeCobras(ii)));
    CobraConfig.ARM_DATA.(fldID).KINEMATICS.Link2_Link_Length = tlMean;
    tlStd = std(TLA); % Take standard deviation of links. Filtering out those greater than 60px
    figure(hf3)
    subplot(9,2,2*ii)
    rlh = refline(0,tlMean,0,'k--');
%     legend(rlh,strcat('Overall Mean=',num2str(tlMean),' (\sigma=',num2str(tlStd),')'));
    legend(rlh,sprintf('Overall Mean=%4.2f \\sigma=%4.2f',tlMean,,tlStd));
end

configSavePath = 'CobraConfig_060914_3.mat';
save(configSavePath,'CobraConfig')
