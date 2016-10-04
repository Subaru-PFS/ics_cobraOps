clear all
close all

% pixel scale
pxSc = 90.5; %um/px

% Find home position logs
homeFiles = ls('home_mId_*pId*.txt');
antihomeFiles = ls('antihome_mId_*pId*.txt');
mId = cellfun(@str2num,(regexp(reshape(homeFiles',1,[]),'(?<=mId_)\d*','match')'));
mpId = cellfun(@str2num,(regexp(reshape(homeFiles',1,[]),'(?<=pId_)\d*','match')'));
pId = (mId-1)*57 + mpId;

% Create home rpt figure
fh1 = figure('name','Home and Antihome Scatters');
fh2 = figure('name','Home Scatters');

% For each log file
for ii=1:length(homeFiles(:,1))

    % Store module and positioner ID in sctuct of LD array
    LD(ii).mId = mId(ii);
    LD(ii).pId = pId(ii);
    
    % Read in raw log data
    LD(ii).rawHome = dlmread(homeFiles(ii,:),',');
    LD(ii).rawAntihome = dlmread(antihomeFiles(ii,:),',');
    
    % Assign raw data to their appropriate fields
    LD(ii).homeP = LD(ii).rawHome(:,1) + LD(ii).rawHome(:,2)*i;
    LD(ii).antihomeP = LD(ii).rawAntihome(:,1) + LD(ii).rawAntihome(:,2)*i;
    LD(ii).homeTheta = LD(ii).rawHome(:,3);
    LD(ii).homePhi = LD(ii).rawHome(:,4);
    LD(ii).antihomeTheta = LD(ii).rawAntihome(:,3);
    LD(ii).antihomePhi = LD(ii).rawAntihome(:,4);
    
    % Check PID7 Phi stage
    pid7phithresh = 8;
    if LD(ii).pId == 7
%         fltrInd = LD(ii).homePhi < pid7phithresh);
        fltrInd = 1:50;
        LD(ii).homePhi = LD(ii).homePhi(fltrInd);
        LD(ii).homeTheta = LD(ii).homeTheta(fltrInd);
        LD(ii).homeP = LD(ii).homeP(fltrInd);
    end
            
    
    % Calculate means and standard deviations
    LD(ii).homeMean = mean(LD(ii).homeP);
    LD(ii).antihomeMean = mean(LD(ii).antihomeP);
    LD(ii).homeStd = std(LD(ii).homeP)/sqrt(2);
    LD(ii).antihomeStd = std(LD(ii).antihomeP)/sqrt(2);
    
    % Create filtered datasets omitting outliers beyond 3sigma
    antiHome_d2mean = abs(LD(ii).antihomeP - LD(ii).antihomeMean);
    antiHomeFltMean = mean(LD(ii).antihomeP(antiHome_d2mean <= 3*LD(ii).antihomeStd));
    antiHomeFltStd = std(LD(ii).antihomeP(antiHome_d2mean <= 3*LD(ii).antihomeStd))/sqrt(2);
    LD(ii).antiHomeFltMean = antiHomeFltMean;
    LD(ii).antiHomeFltStd = antiHomeFltStd;
    LD(ii).antiHome_d2mean = antiHome_d2mean;
    home_d2mean = abs(LD(ii).homeP - LD(ii).homeMean);
    homeFltMean = mean(LD(ii).homeP(home_d2mean <= 3*LD(ii).homeStd));
    homeFltStd = std(LD(ii).homeP(home_d2mean <= 3*LD(ii).homeStd))/sqrt(2);
    LD(ii).homeFltMean = homeFltMean;
    LD(ii).homeFltStd = homeFltStd;
    LD(ii).home_d2mean = home_d2mean;
        
    
    % Create plot of positions and statistics
    figure(fh1)
%     ph1 = plot([LD(ii).homeP,LD(ii).antihomeP],'x');
    ph1(1) = plot(LD(ii).homeP,'bx');
    hold on;
    ph1(2) = plot(LD(ii).antihomeP,'gx');
    axis equal;
    ph2 = cmplx(@circle,LD(ii).homeMean,3*LD(ii).homeStd,'b');
    th2 = cmplx(@text,LD(ii).homeMean,sprintf('3\\sigma circle (\\sigma=%4.2fum)',LD(ii).homeStd*pxSc),'color','b');
    ph3 = cmplx(@circle,LD(ii).antihomeMean,3*LD(ii).antihomeStd,'g');
    th3 = cmplx(@text,LD(ii).antihomeMean,sprintf('3\\sigma circle (\\sigma=%4.2fum)',LD(ii).antihomeStd*pxSc),'color','g');
    
    ph4 = cmplx(@circle,LD(ii).antiHomeFltMean,3*LD(ii).antiHomeFltStd,'k');
    th4 = cmplx(@text,LD(ii).antiHomeFltMean,sprintf('3\\sigma circle (\\sigma=%4.2fum)',LD(ii).antiHomeFltStd*pxSc),'color','k');
    ph5 = cmplx(@circle,LD(ii).homeFltMean,3*LD(ii).homeFltStd,'k');
    th5 = cmplx(@text,LD(ii).homeFltMean,sprintf('3\\sigma circle (\\sigma=%4.2fum)',LD(ii).homeFltStd*pxSc),'color','k');
    
    figure(fh2)
    sph = subplot(3,3,ii);
    cpp = copyobj([ph1(1),ph2,ph5],sph);
    axis equal
    if LD(ii).homeFltStd*pxSc < 10
        RorG = 'g';
    else
        RorG = 'r';
    end
    title(sprintf('pid%d',ii),'color',RorG)
    legend(sprintf('Filtered 3\\sigma circle (\\sigma=%4.2fum)',LD(ii).homeFltStd*pxSc),sprintf('3\\sigma circle (\\sigma=%4.2fum)',LD(ii).homeStd*pxSc))
    
    
end

CW = lines(length(LD));

fh3 = figure;
subplot(2,1,1)
for ii=1:length(LD)
    plot(LD(ii).home_d2mean,'color',CW(ii,:))
    hold on
end
% plot([LD.home_d2mean])
legend(cellfun(@num2str,{LD.pId},'UniformOutput',0))
title('Distance from Mean Home')
ylabel('px')
xlabel('trial')
subplot(2,1,2)
for ii=1:length(LD)
    plot(LD(ii).antiHome_d2mean,'color',CW(ii,:))
    hold on
end
% plot([LD.antiHome_d2mean])
title('Distance from Mean Antihome')
ylabel('px')
xlabel('trial')

fh4 = figure;
subplot(2,1,1)
for ii=1:length(LD)
    plot(LD(ii).homeTheta,'color',CW(ii,:))
    hold on
end
% plot([LD.homeTheta])
title('Home Theta')
ylabel('theta')
xlabel('trial')
legend(cellfun(@num2str,{LD.pId},'UniformOutput',0))
subplot(2,1,2)
for ii=1:length(LD)
    plot(LD(ii).homePhi,'color',CW(ii,:))
    hold on
end
% plot([LD.homePhi])
title('Home Phi')
ylabel('phi')
xlabel('trial')

fh5 = figure;
subplot(2,1,1)
for ii=1:length(LD)
    plot(LD(ii).antihomeTheta,'color',CW(ii,:))
    hold on
end
% plot([LD.antihomeTheta])
title('antihome Theta')
ylabel('theta')
xlabel('trial')
legend(cellfun(@num2str,{LD.pId},'UniformOutput',0))
subplot(2,1,2)
plot([LD.antihomePhi])
title('antihome Phi')
ylabel('phi')
xlabel('trial')
















