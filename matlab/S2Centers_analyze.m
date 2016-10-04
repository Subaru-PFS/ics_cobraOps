% $$$ % S1CentersGroupAnlys
% $$$ 
% $$$ clear all
% $$$ close all
% $$$ 
% $$$ c=1;
% $$$ imgsFiles = {};
% $$$ dataTitle = {};
% $$$ 
% $$$ % %3-28-14 S2 Centers
% $$$ % imgsFiles{c} = '..\..\TEST_RESULTS\Metrology\centerms2_033114_1.mat';
% $$$ % dataTitle{c} = '3/28/14-1';
% $$$ % c=c+1;
% $$$ 
% $$$ metDir = '..\..\TEST_RESULTS\Metrology\';
% $$$ metDirFiles = dir2cell(metDir);
% $$$ phiFileInd = strmatch('phiMtrlgy',metDirFiles);
% $$$ phiFiles = {metDirFiles{phiFileInd}};
% $$$ imgsFiles = [imgsFiles fullfile(metDir,phiFiles)];
% $$$ dataTitle = [dataTitle phiFiles];



%% EXECUTION
%close all
clear all

MAXLINK = 60; % pixels?
load phiMtrlgySummary_051414 data dataTitle imgsFiles
data = cell2mat(data);

for jj = 1:length(data)
  %% filter outliers
  tempL1s = data(jj).Link1s;
  tempL2s = data(jj).Link2s;
  tempL1s(tempL1s > MAXLINK) = NaN;
  tempL2s(tempL2s > MAXLINK) = NaN;
  
  
  %% make new fields for easy concatenation
  data(jj).L1Means = nanmean(tempL1s);
  data(jj).L2Means = nanmean(tempL2s);
  tdate = datenum(regexp(imgsFiles{jj},'\d+_\d+_\d+','match'),'mm_dd_yy');
  data(jj).dates = repmat(tdate, size(tempL1s,1), 1);
  data(jj).date  = tdate;
  
  clear tempL1s, tempL2s;
end

Link1s  = vertcat(data.Link1s);
Link2s  = vertcat(data.Link2s);
L1Means = vertcat(data.L1Means);
L2Means = vertcat(data.L2Means);
dates   = vertcat(data.dates); % matches LinkNs
date    = vertcat(data.date);  % matches LNMeans

fh1 = figure(1);%'name','Link1s vs days');
fh2 = figure(2);%'name','Link2s vs days');

for jj = 1:size(data(jj).Link1s,2)
  figure(fh1)
  subplot(2,3,jj)
  plot(dates - min(dates), Link1s(:,jj) , 'ro'); hold on;
  plot(date - min(date)  , L1Means(:,jj), 'bs'); hold off;
  xlabel('days'); ylabel('length [pix]');
  title(sprintf('pId%d',jj));
  ylim([35 50]);

  figure(fh2)
  subplot(2,3,jj)
  plot(dates - min(dates), Link2s(:,jj) , 'bo'); hold on;
  plot(date - min(date)  , L2Means(:,jj), 'rs'); hold off;
  xlabel('days'); ylabel('length [pix]');
  title(sprintf('pId%d',jj));
  ylim([35 50]);
end

return
figure(3)


hf3 = existingFigure('Link Lengths');
for ii=1:length(Link1s(1,:))
    TLA = Link1s(:,ii);% Three letter acronym meaning temporary link array
    tlMean = mean(TLA(TLA<60)); % Take mean of links. Filtering out those greater than 60px
    tlStd = std(TLA(TLA<60)); % Take standard deviation of links. Filtering out those greater than 60px
    figure(hf3)
    subplot(5,2,2*(ii-1)+1)
    rlh = refline(0,tlMean,0,'k--');
    legend(rlh,strcat('Overall Mean=',num2str(tlMean),' (\sigma=',num2str(tlStd),')'));
end

for ii=1:length(Link2s(1,:))
    TLA = Link2s(:,ii);% Three letter acronym meaning temporary link array
    tlMean = mean(TLA(TLA<60)); % Take mean of links. Filtering out those greater than 60px
    tlStd = std(TLA(TLA<60)); % Take standard deviation of links. Filtering out those greater than 60px
    figure(hf3)
    subplot(5,2,2*ii)
    rlh = refline(0,tlMean,0,'k--');
    legend(rlh,strcat('Overall Mean=',num2str(tlMean),' (\sigma=',num2str(tlStd),')'));
end
    

% for ii=1:length
% % Calculate overall mean center
% cMean = mean(cMeans);
% 
% % Calculate std deviation of means
% cStd = std(cMeans);
% 
% figure(1)
% ch = cmplx(@circle,cMean,3*cStd,'k');
% appendLegend(ch2,strcat('All Points -- 3\sigma circle (\sigma=',num2str(cStd),')'))

    
% OverallMean = 
CW = hsv(length(dataTitle));
figure
L1N = [];
for ii=1:length(dataTitle)-1
    load(dataTitle{ii})
    for jj=1:length(data{ii}.dnames)
%         keyboard
        thisImg = eval(data{ii}.dnames{jj});
        thisN = [];
        for ll=1:5
            fldID = strcat('pId',num2str(ll));
            thisN(ll) = length(thisImg.(fldID).CCDpos);
        end
        L1N = vertcat(L1N,thisN);
        
        
        plotFibers(thisImg,'x','color',CW(ii,:))
        hold on
    end
end

            
            
        
