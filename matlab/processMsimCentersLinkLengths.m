% processMsimMaps
clear all
close all

%% INPUTS
mid = 1;
% Provide the config file with which the data was taken with
cfgfiles = dir2cell('*.xml');
if length(cfgfiles) == 1
  cli_answer = input(sprintf('Use %s? [Y|n] ',cfgfiles{1}),'s');
  if strcmp(cli_answer,'n')
    baseCfg = loadCfgXml;
  else
    baseCfg = loadCfgXml('.',cfgfiles{1});
  end
else
  baseCfg = loadCfgXml;
end

  
%% EXECUTION (Don't edit below here)
% Get theta dots.
thetaFw = makeMsimCenters(1,1,'thetaFwMap');

% Get phi dots.
phiFw = makeMsimCenters(2,1,'PhiFwMap');

% Start new config
newConfig = baseCfg;
data = [];
c=1;
for pid = [thetaFw.pId]
    pidname = sprintf('pid%d', pid);
    
   for kk = 1:sum(regexpcmp(fieldnames(baseCfg.ARM_DATA),'ARM_DATA_*')) % Count all struct that start with the name ARM_DATA_
        ARM_DATA_n = sprintf('ARM_DATA_%d', kk);
        if(pid ==  str2num(baseCfg.ARM_DATA.(ARM_DATA_n).DATA_HEADER.Positioner_Id.Text))
            armid = ARM_DATA_n;
            break;
        end
   end
    
    thetaPositions = [];
    thetaPositions =  thetaFw(c).rawData;
    thetaPositions = thetaPositions(:,1:2);
   
    thetaCfit0 = circfit(thetaPositions(:,1),thetaPositions(:,2));
    
    figure(10)
    axis([0 2000 0 2000]);
    daspect([1,1,1]);
    hold on
    plot(thetaCfit0.xc, thetaCfit0.yc,'rx');
    plot(thetaPositions(:,1),thetaPositions(:,2), 'bo');
   
    thetaCentroids = thetaPositions * [1;1i];
    thetaCentroidLengths = thetaCentroids - thetaCfit0.c;
    thetaRadii = abs(thetaCentroidLengths);
    thetaAngles = angle(thetaCentroidLengths);
    thetaUangles = unwrap(thetaAngles); % use a continuous set of angles for calcs
    
    % filter out hard stops.  need to deal with the fact that diff
    % reduces the length of the vector by 1.
    thetaOk = (abs(diff(thetaUangles)) > .01);
    thetaOk = thetaOk & [0; thetaOk(1:end-1)];
    % refit the data    
    cfit1 = circfit(thetaCentroids(thetaOk));
    thetaCentroidLengths = thetaCentroids - cfit1.c;
    thetaRadii = abs(thetaCentroidLengths);
    thetaAngles = angle(thetaCentroidLengths);

    figure(11) 
    plot(thetaAngles/(2*pi),thetaRadii, 'bx');
    hold on;
    plot(thetaAngles(thetaOk)/(2*pi), thetaRadii(thetaOk),'go');%,'MarkerFace','g');
    hold off;
    xlabel('angle/\tau');
    ylabel('radius [pix]');
    title(sprintf('THETA PID %d',pid));
    figure(12)
    hist(thetaRadii(thetaOk), 50)

    data.(pidname).thetaRadius = thetaRadii(thetaOk);
    data.(pidname).thetaCenter  = cfit1.c;
    disp('theta');
    %keyboard;
    
    
    
    
    phiPositions =  phiFw(c).rawData;
    phiPositions = phiPositions(:,1:2);
   
    phiCfit0 = circfit(phiPositions(:,1),phiPositions(:,2));
    
    figure(20)
    axis([0 2000 0 2000]);
    daspect([1,1,1]);
    hold on
    plot(phiCfit0.xc, phiCfit0.yc,'rx');
    plot(phiPositions(:,1),phiPositions(:,2), 'bo');
   
    phiCentroids = phiPositions * [1;1i];
    phiCentroidLengths = phiCentroids - phiCfit0.c;
    phiRadii = abs(phiCentroidLengths);
    phiAngles = angle(phiCentroidLengths);
    phiUangles = unwrap(phiAngles); % use a continuous set of angles for calcs
    
    % filter out hard stops.  need to deal with the fact that diff
    % reduces the length of the vector by 1.
    phiOk = (abs(diff(phiUangles)) > .01);
    phiOk = phiOk & [0; phiOk(1:end-1)];
    % refit the data    
    cfit1 = circfit(phiCentroids(phiOk));
    phiCentroidLengths = phiCentroids - cfit1.c;
    phiRadii = abs(phiCentroidLengths);
    phiAngles = angle(phiCentroidLengths);

    figure(21) 
    plot(phiAngles/(2*pi),phiRadii, 'bx');
    hold on;
    plot(phiAngles(phiOk)/(2*pi), phiRadii(phiOk),'go');%,'MarkerFace','g');
    hold off;
    xlabel('angle/\tau');
    ylabel('radius [pix]');
    title(sprintf('PHI PID %d',pid));
    figure(22)
    hist(phiRadii(phiOk), 50)

    data.(pidname).phiRadius = phiRadii(phiOk);
    data.(pidname).phiCenter  = cfit1.c;
    disp('phi'); 
    
    keyboard;
    
    
    data.(pidname).phiRadius = phiRadii(phiOk);
    data.(pidname).phiCenter  = cfit1.c;
    %% Now the calculate the link lengths:
    
    

       c=c+1;
end

% Take the last arm data structure and store it as filler for dummy arms
ARM_DATA_FILLER = newConfig.ARM_DATA.(armid);

% Create dummy arm data structures for missing cobras
for pid = 1:max([thetaFw.pId])
    if isempty(find([thetaFw.pId]==pid))
        armid = sprintf('ARM_DATA_%d',pid);
        newConfig.ARM_DATA.(armid) = ARM_DATA_FILLER;
    end
end

% % Save new XML with adjusted maps
% [xmlfile, xmlfilepath] = uiputfile('*.xml','Save new CobraConfig XML file with new motor maps');
% cobraCfg2xml(newConfig,fullfile(xmlfilepath,xmlfile));
%