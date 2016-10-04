%TARGET GENERATOR
clear all
close all

%% INPUTS
% Get config xml
CobraConfig = loadCfgXml;
% CobraConfig = loadCfgXml('C:\Users\sage\Desktop\Dropbox\PFS_EM\TEST_RESULTS\Metrology','LinkLengths_Centers_062614_2.xml');
%CobraConfig2 = loadCfgXml;
% Specify which cobras to create targets for
moduleID = 1;
%activeCobras = [1,5,9,13,17,21,25,29,33,37,41,45,49,53]; % Same purpose as movingCobras but more intuitive to use pId
%activeCobras = [9,15,31,37,45]; % Same purpose as movingCobras but more intuitive to use pId
%activeCobras = [3, 7, 11, 15, 19, 23, 27, 31, 35, 39, 43, 47, 51];
activeCobras = [1, 3, 5, 7, 9, 11, 13, 15, 17, 19, 21, 23, 25, 27, 29, 31, 33, 35, 37, 39, 41, 43, 45, 47, 49, 51, 53];

% Specify number of random targets to generate
numtrg = 200;

reachAround = 0.05; % Specify reach factor. Eg: 0.05 means that 5% of 
% reachable radial distance is excluded for target generation so that there
% is some padding on ideal keepout zones.

LF = [10*ones(1,numtrg-1) 48]; % stupid hack to remove \n after last target.

%% EXECUTION
fh1 = figure;
hold on;
kl = 1;

for ii = activeCobras
  cnfgID = sprintf('ARM_DATA_%d',ii);
  
  tfile = fopen(sprintf('TargetList_mId_%d_pId_%d.txt',moduleID,ii),'w');
  
  center = getARMval(CobraConfig,ii, moduleID,'Global_base_pos_x') + ...
           getARMval(CobraConfig,ii, moduleID, 'Global_base_pos_y')*i;
  L1 = getARMval(CobraConfig,ii, moduleID,'Link1_Link_Length');
  L2 =  getARMval(CobraConfig,ii, moduleID,'Link2_Link_Length');
  phiMax = getARMval(CobraConfig,ii, moduleID,'phiMax');
  phiMin = getARMval(CobraConfig,ii, moduleID,'phiMin');

  if isempty(phiMax)
    phiMax = 160;
  end
  if isempty(phiMin)
    phiMin = 5;
  end
  % Inner keepout radius
  Riko = sqrt(L1^2+L2^2-2*L1*L2*cos(phiMin*pi/180));
  % Max reach
  reach = sqrt(L1^2+L2^2-2*L1*L2*cos(phiMax*pi/180));
  
  if(kl==1)
    cmplx(@text,center + 20 + 20*1i,num2str(ii));
    plot(center,'r.')
    cmplx(@circle,center,reach,'r');
    
  elseif(kl==2)
    cmplx(@text,center + 20 - 20*1i,num2str(ii)); 
    plot(center,'g.')
    cmplx(@circle,center,reach,'g');
  end

  %% generate target locations
  gamma = rand(1,numtrg)*2*pi;
  pad = reachAround*reach;
  bad  = 1:numtrg; % "bad" radii are outside the range 
                   % Riko + pad : reach - pad
  while ~isempty(bad)
    Rtrg(bad) = sqrt(rand(size(bad))) * (reach - pad);
    bad = find(Rtrg < (Riko + pad));
  end
  
  target = center + Rtrg .* exp(i*gamma);
  plot(target,'r.');
  % label targets with integers
  cmplx(@text,target,strsplit(num2str(1:length(target)),' '));
  drawnow;
  
  targetoutput = [simple(target.')' ; LF];
  
  fprintf(tfile,'%4.2f,%4.2f%c',targetoutput);
  
% $$$      for ii=1:numtrg
% $$$         % Pick a random angle between 0 and 2pi
% $$$         gamma = rand*2*pi;
% $$$         phi = phiMin*pi/180+ rand * (pi-(phiMax+phiMin)*pi/180);
% $$$         % Pick a random radii reachable by arm
% $$$         pad = reachAround*reach;
% $$$         Rtrg = Riko + pad + rand*(reach-Riko-2*pad);
% $$$ 	% Rtrg = sqrt(rand * (reach - pad));
% $$$ 
% $$$ 	
% $$$         % Targets equally spread out over the angle:
% $$$ %         tar = L1 * exp(1i * gamma) + L2 * exp(1i * 1);
% $$$ %         target = center + tar; 
% $$$ %         
% $$$         % Calculate x and y components of target from cobra center
% $$$         cxT = Rtrg*cos(gamma);
% $$$         cyT = Rtrg*sin(gamma);
% $$$         
% $$$         % Add cobra center to get target x,y in fiducial frame
% $$$         target = center + cxT + cyT*i;
% $$$         
% $$$         plot(target,'ro','MarkerFaceColor','red','MarkerSize',1.5);
% $$$         cmplx(@text,target,num2str(ii));
% $$$         
% $$$         if ii==numtrg
% $$$             fprintf(tfile,'%4.2f,%4.2f',real(target),imag(target));
% $$$         else
% $$$             fprintf(tfile,'%4.2f,%4.2f\n',real(target),imag(target));
% $$$         end
% $$$     end
  
  %     target = center + L1 + L2*i;
  %     fprintf(tfile,'%4.2f,%4.2f',real(target),imag(target));
  
  fclose(tfile);

end

axis equal
hold off;
